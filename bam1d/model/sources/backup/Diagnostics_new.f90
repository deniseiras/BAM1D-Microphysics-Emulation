!
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:37 $
!  $Revision: 1.22 $
!
MODULE Diagnostics
  USE Constants, ONLY :          &
       r8,i8,r4,i4,&
       grav   , &
       tov  

  USE Options, ONLY: &
       nfprt,            &
       nferr,            &
       ifprt,            &
       yrl  , &
       monl , &
       cthl
 
  USE IOLowLevel, ONLY: &
       WriteField   , &
       WriteProgHead, &
       WriteDir     , &
       WriteDire 

  USE InputOutput, ONLY: &
       transp,           &
       sclout,           &
       cnvray
 
  USE Utils, ONLY: &
       tmstmp2

  USE Constants, ONLY :   &
       r8, i8, &
       ndavl, ndrq, ncdg, jxavl, jxcdg, numx, grav,numx


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: InitDiagnostics
  PUBLIC :: pwater
  PUBLIC :: rsdiag
  PUBLIC :: updia
  PUBLIC :: wridia
  PUBLIC :: wrprog

  PUBLIC :: lgaus
  PUBLIC :: combf
  PUBLIC :: reqdg
  PUBLIC :: itcf
  PUBLIC :: nucf
  PUBLIC :: lvcf
  PUBLIC :: nurq
  PUBLIC :: iavrq
  PUBLIC :: itavl
  PUBLIC :: nuavl
  PUBLIC :: lvavl
  PUBLIC :: dodia
  PUBLIC :: inavl
  PUBLIC :: ixavl
  PUBLIC :: iclcd
  PUBLIC :: incf
  PUBLIC :: ixcf
  PUBLIC :: kravl
  PUBLIC :: jrcf
  PUBLIC :: krcf
  PUBLIC :: icf
  PUBLIC :: mxavl
  PUBLIC :: ngaus
  PUBLIC :: mgaus
  PUBLIC :: gaus  
  PUBLIC :: gaus_in
  PUBLIC :: lvrq

  INTEGER, ALLOCATABLE :: SMap(:)
  INTEGER, ALLOCATABLE :: AMap(:)
  INTEGER :: sMapSize
  INTEGER :: aMapSize


  INTEGER, PARAMETER :: ndp=100


  REAL(KIND=r8),    PARAMETER :: eeemin=1.0e-35_r8
  REAL(KIND=r8),    PARAMETER :: eeemax=1.0e35_r8
  REAL(KIND=r8),    PARAMETER :: undef =1.0e53_r8
  INTEGER, PARAMETER :: ngbme = 13
  INTEGER, PARAMETER :: igbme(ngbme)= &
                               (/9,10,16,18,19,17,21,23,24,25,28,29,30/)
  CHARACTER(LEN=4), PARAMETER :: ivar(2) = &
                                          (/'GAUS', 'GAUS'/)
  CHARACTER(LEN=256) :: fname
  INTEGER :: ierr
 
!--------------------------------------------------------------------------
!         locations for available diagnostics in this subroutine
!--------------------------------------------------------------------------     
  ! Available Diagnostics Indexes
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmpsfc =  1 ! time mean surface pressure
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmdivg =  2 ! time mean divergence (subroutine accpf)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmvort =  3 ! time mean vorticity (subroutine accpf)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmsphu =  4 ! time mean specific humidity (subroutine accpf)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmtvir =  5 ! time mean virtual temperature (subroutine accpf)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmtsfc =  6 ! time mean surface temperature
  INTEGER, PUBLIC, PARAMETER :: nDiag_omegav =  7 ! omega
  INTEGER, PUBLIC, PARAMETER :: nDiag_sigdot =  8 ! sigma dot
  INTEGER, PUBLIC, PARAMETER :: nDiag_toprec =  9 ! total precipiation
  INTEGER, PUBLIC, PARAMETER :: nDiag_cvprec = 10 ! convective precipitation
  INTEGER, PUBLIC, PARAMETER :: nDiag_lsprec = 11 ! large scale precipitation
  INTEGER, PUBLIC, PARAMETER :: nDiag_snowfl = 12 ! snowfall
  INTEGER, PUBLIC, PARAMETER :: nDiag_runoff = 13 ! runoff
  INTEGER, PUBLIC, PARAMETER :: nDiag_pwater = 14 ! precipitable water
  INTEGER, PUBLIC, PARAMETER :: nDiag_intlos = 15 ! interception loss
  INTEGER, PUBLIC, PARAMETER :: nDiag_sheatf = 16 ! sensible heat flux
  INTEGER, PUBLIC, PARAMETER :: nDiag_lheatf = 17 ! latent heat flux
  INTEGER, PUBLIC, PARAMETER :: nDiag_ustres = 18 ! surface zonal stress
  INTEGER, PUBLIC, PARAMETER :: nDiag_vstres = 19 ! surface meridional stress
  INTEGER, PUBLIC, PARAMETER :: nDiag_cloudc = 20 ! cloud cover
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwdbot = 21 ! longwave downward at bottom
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwubot = 22 ! longwave upward at bottom
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwutop = 23 ! longwave upward at top
  INTEGER, PUBLIC, PARAMETER :: nDiag_swdtop = 24 ! shortwave downward at top
  INTEGER, PUBLIC, PARAMETER :: nDiag_swdbot = 25 ! shortwave downward at ground
  INTEGER, PUBLIC, PARAMETER :: nDiag_swubot = 26 ! shortwave upward at bottom
  INTEGER, PUBLIC, PARAMETER :: nDiag_swutop = 27 ! shortwave upward at top
  INTEGER, PUBLIC, PARAMETER :: nDiag_swabea = 28 ! shortwave absorbed by the earth/atmosphere
  INTEGER, PUBLIC, PARAMETER :: nDiag_swabgr = 29 ! shortwave absorbed by the ground
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwnetb = 30 ! net longwave at bottom
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwheat = 31 ! longwave heating
  INTEGER, PUBLIC, PARAMETER :: nDiag_swheat = 32 ! shortwave heating
  INTEGER, PUBLIC, PARAMETER :: nDiag_clheat = 33 ! convective latent heating
  INTEGER, PUBLIC, PARAMETER :: nDiag_cmchan = 34 ! convective moisture change
  INTEGER, PUBLIC, PARAMETER :: nDiag_lslhea = 35 ! large scale latent heating
  INTEGER, PUBLIC, PARAMETER :: nDiag_lsmcha = 36 ! large scale moisture change
  INTEGER, PUBLIC, PARAMETER :: nDiag_sclhea = 37 ! shallow convective latent heating
  INTEGER, PUBLIC, PARAMETER :: nDiag_scmcha = 38 ! shallow convective moisture change
  INTEGER, PUBLIC, PARAMETER :: nDiag_vdheat = 39 ! vertical diffusion heating
  INTEGER, PUBLIC, PARAMETER :: nDiag_vdmois = 40 ! vertical diffusion moistening
  INTEGER, PUBLIC, PARAMETER :: nDiag_vduzon = 41 ! vertical diffusion zonal momentum change
  INTEGER, PUBLIC, PARAMETER :: nDiag_vdvmer = 42 ! vertical diffusion meridional momentum change
  INTEGER, PUBLIC, PARAMETER :: nDiag_txgwds = 43 ! gravity wave drag surface zonal stress
  INTEGER, PUBLIC, PARAMETER :: nDiag_tygwds = 44 ! gravity wave drag surface meridional stress
  INTEGER, PUBLIC, PARAMETER :: nDiag_gwduzc = 45 ! gravity wave drag zonal momentum change
  INTEGER, PUBLIC, PARAMETER :: nDiag_gwdvmc = 46 ! gravity wave drag meridional momentum change
  INTEGER, PUBLIC, PARAMETER :: nDiag_hhedif = 47 ! horizontal heating diffusion
  INTEGER, PUBLIC, PARAMETER :: nDiag_hmodif = 48 ! horizontal moisture diffusion
  INTEGER, PUBLIC, PARAMETER :: nDiag_hdidif = 49 ! horizontal divergence diffusion
  INTEGER, PUBLIC, PARAMETER :: nDiag_hvodif = 50 ! horizontal vorticity diffusion
  INTEGER, PUBLIC, PARAMETER :: nDiag_divgxq = 51 ! divergence * specific humidity
  INTEGER, PUBLIC, PARAMETER :: nDiag_vmoadv = 52 ! vertical moisture advection
  INTEGER, PUBLIC, PARAMETER :: nDiag_hmofcv = 53 ! horizontal moisture flux convergence (???)
  INTEGER, PUBLIC, PARAMETER :: nDiag_vimfcv = 54 ! vertically integrated moisture flux convergence (???)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmlnps = 55 ! time mean log surface pressure (subroutine accpf)
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwdbtc = 56 ! longwave downward at bottom (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwutpc = 57 ! longwave upward at top (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_swdbtc = 58 ! shortwave downward at ground (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_swubtc = 59 ! shortwave upward at bottom (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_swutpc = 60 ! shortwave upward at top (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_swaeac = 61 ! shortwave absorbed by the earth/atmosphere (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_swabgc = 62 ! shortwave absorbed by the ground (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwnbtc = 63 ! net longwave at bottom (clear)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tmtdps = 64 ! time mean deep soil temperature
  INTEGER, PUBLIC, PARAMETER :: nDiag_tgfccv = 65 ! ground/surface cover temperature
  INTEGER, PUBLIC, PARAMETER :: nDiag_tcanop = 66 ! canopy temperature
  INTEGER, PUBLIC, PARAMETER :: nDiag_tcairs = 67 ! temperature of canopy air space
  INTEGER, PUBLIC, PARAMETER :: nDiag_ecairs = 68 ! vapor pressure of canopy air space
  INTEGER, PUBLIC, PARAMETER :: nDiag_bsolht = 69 ! bare soil latent heat
  INTEGER, PUBLIC, PARAMETER :: nDiag_nshcrm = 70 ! negative specific humidity correction moisture source
  INTEGER, PUBLIC, PARAMETER :: nDiag_ozonmr = 71 ! ozone mass mixing ratio (g/g)
  INTEGER, PUBLIC, PARAMETER :: nDiag_vdtclc = 72 ! vertical dist total cloud cover
  INTEGER, PUBLIC, PARAMETER :: nDiag_invcld = 73 ! inversion cloud
  INTEGER, PUBLIC, PARAMETER :: nDiag_ssatcl = 74 ! supersaturation cloud
  INTEGER, PUBLIC, PARAMETER :: nDiag_cnvcld = 75 ! convective cloud
  INTEGER, PUBLIC, PARAMETER :: nDiag_shcvcl = 76 ! shallow convective cloud
  INTEGER, PUBLIC, PARAMETER :: nDiag_clliwp = 77 ! cloud liquid water path
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwcemi = 78 ! longwave cloud emissivity
  INTEGER, PUBLIC, PARAMETER :: nDiag_sclopd = 79 ! shortwave cloud optical depth
  INTEGER, PUBLIC, PARAMETER :: nDiag_mofres = 80 ! momentum flux resistance
  INTEGER, PUBLIC, PARAMETER :: nDiag_casrrs = 81 ! canopy air spc to ref. lvl resistance
  INTEGER, PUBLIC, PARAMETER :: nDiag_cascrs = 82 ! canopy air spc to canopy resistance
  INTEGER, PUBLIC, PARAMETER :: nDiag_casgrs = 83 ! canopy air spc to ground resistance
  INTEGER, PUBLIC, PARAMETER :: nDiag_gcovrs = 84 ! ground cover resistance
  INTEGER, PUBLIC, PARAMETER :: nDiag_bssfrs = 85 ! bare soil surface resistance
  INTEGER, PUBLIC, PARAMETER :: nDiag_homtvu = 86 ! Horizontal Momentum Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_vzmtwu = 87 ! Vertical Zonal Momentum Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_vmmtwv = 88 ! Vertical Meridional Momentum Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_mshtvt = 89 ! Meridional Sensible Heat Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_zshtut = 90 ! Zonal Sensible Heat Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_vshtwt = 91 ! Vertical Sensible Heat Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_mshtuq = 92 ! Meridional Specific Humidity Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_zshtuq = 93 ! Zonal Specific Humidity Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_vshtwq = 94 ! Vertical Specific Humidity Transport
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwhtcl = 95 ! longwave heating (clear sky)
  INTEGER, PUBLIC, PARAMETER :: nDiag_swhtcl = 96 ! shortwave heating (clear sky)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tep02m = 97 ! time mean temp at 2-m from sfc layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_mxr02m = 98 ! time mean es humid at 2-m from sfc layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_spw02m = 99 ! Speed wind at 2-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_tep10m = 100 ! Temperature at 10-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_mxr10m = 101 ! specifc humidity at 10-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_spw10m = 102 ! Speed wind at 10-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_viozoc = 103 ! Vertically Integrated Ozone Content (Dobson units)
  INTEGER, PUBLIC, PARAMETER :: nDiag_dewptt = 104 ! Dew Point Temperature K
  INTEGER, PUBLIC, PARAMETER :: nDiag_zwn10m = 105 ! Zonal Wind at 10-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_mwn10m = 106 ! Meridional wind at 10-m from surface layer
  INTEGER, PUBLIC, PARAMETER :: nDiag_biomau = 107 ! Total biomass in the upper canopy (kg_C m-2)
  INTEGER, PUBLIC, PARAMETER :: nDiag_biomal = 108 ! Total biomass in the lower canopy (kg_C m-2)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tlaiup = 109 ! Total leaf area index for the upper canopy
  INTEGER, PUBLIC, PARAMETER :: nDiag_tlailw = 110 ! Total leaf area index for the lower canopy
  INTEGER, PUBLIC, PARAMETER :: nDiag_tstnsp = 111 ! Total storage of N in soil profile (kg_N m-2) 
  INTEGER, PUBLIC, PARAMETER :: nDiag_wsttot = 112 ! Total amount of water stored in snow, soil, puddels, and on vegetation (kg_h2o)
                                                   ! fraction of root in soil layer 
  INTEGER, PUBLIC, PARAMETER :: nDiag_lidecf = 113 ! litter decomposition factor                  (dimensionless)
  INTEGER, PUBLIC, PARAMETER :: nDiag_somdfa = 114 ! soil organic matter decomposition factor        (dimensionless)
  INTEGER, PUBLIC, PARAMETER :: nDiag_facuca = 115 ! frac overall area cover by upper canopy 
  INTEGER, PUBLIC, PARAMETER :: nDiag_fsfclc = 116 ! frac snowfree area cover by lower canopy
  INTEGER, PUBLIC, PARAMETER :: nDiag_frsnow = 117 ! fractional snow cover
  INTEGER, PUBLIC, PARAMETER :: nDiag_insnpp = 118 ! instantaneous npp (mol-CO2 / m-2 / second)
  INTEGER, PUBLIC, PARAMETER :: nDiag_insnee = 119 ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
  INTEGER, PUBLIC, PARAMETER :: nDiag_grbdy0 = 120 ! annual total growing degree days for current year > 0C
  INTEGER, PUBLIC, PARAMETER :: nDiag_grbdy5 = 121 ! annual total growing degree days for current year > 5C  
  INTEGER, PUBLIC, PARAMETER :: nDiag_avet2m = 122 ! monthly average 2-m surface-air temperature 
  INTEGER, PUBLIC, PARAMETER :: nDiag_monnpp = 123 ! monthly total npp for ecosystem (kg-C/m**2/month)
  INTEGER, PUBLIC, PARAMETER :: nDiag_monnee = 124 ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
  INTEGER, PUBLIC, PARAMETER :: nDiag_yeanpp = 125 ! annual total npp for ecosystem (kg-c/m**2/yr)
  INTEGER, PUBLIC, PARAMETER :: nDiag_yeanee = 126 ! annual total NEE for ecosystem (kg-C/m**2/yr) 
  INTEGER, PUBLIC, PARAMETER :: nDiag_upclai = 127 ! upper canopy single-sided leaf area index (area leaf/area veg)
  INTEGER, PUBLIC, PARAMETER :: nDiag_lwclai = 128 ! lower canopy single-sided leaf area index (area leaf/area veg)
  INTEGER, PUBLIC, PARAMETER :: nDiag_pblstr = 129 ! surface friction velocity
  INTEGER, PUBLIC, PARAMETER :: nDiag_hghpbl = 130 ! planetary boundary layer height
  INTEGER, PUBLIC, PARAMETER :: nDiag_khdpbl = 131 ! diffusion coefficient for heat
  INTEGER, PUBLIC, PARAMETER :: nDiag_kmdpbl = 132 ! diffusion coefficient for momentum
  INTEGER, PUBLIC, PARAMETER :: nDiag_ricpbl = 133 ! bulk Richardson no. from level to ref lev

  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts01 = 134 ! pft tropical broadleaf evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts02 = 135 ! pft tropical broadleaf drought-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts03 = 136 ! pft warm-temperate broadleaf evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts04 = 137 ! pft temperate conifer evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts05 = 138 ! pft temperate broadleaf cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts06 = 139 ! pft boreal conifer evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts07 = 140 ! pft boreal broadleaf cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts08 = 141 ! pft boreal conifer cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts09 = 142 ! pft evergreen shrubs
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts10 = 143 ! pft cold-deciduous shrubs
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts11 = 144 ! pft warm (c4) grasses
  INTEGER, PUBLIC, PARAMETER :: nDiag_pfts12 = 145 ! pft cool (c3) grasses

  INTEGER, PUBLIC, PARAMETER :: nDiag_biol01 = 146 ! cbiol tropical broadleaf evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol02 = 147 ! cbiol tropical broadleaf drought-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol03 = 148 ! cbiol warm-temperate broadleaf evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol04 = 149 ! cbiol temperate conifer evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol05 = 150 ! cbiol temperate broadleaf cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol06 = 151 ! cbiol boreal conifer evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol07 = 152 ! cbiol boreal broadleaf cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol08 = 153 ! cbiol boreal conifer cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol09 = 154 ! cbiol evergreen shrubs
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol10 = 155 ! cbiol cold-deciduous shrubs
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol11 = 156 ! cbiol warm (c4) grasses
  INTEGER, PUBLIC, PARAMETER :: nDiag_biol12 = 157 ! cbiol cool (c3) grasses

  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp01 = 158 ! ynpp tropical broadleaf evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp02 = 159 ! ynpp tropical broadleaf drought-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp03 = 160 ! ynpp warm-temperate broadleaf evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp04 = 161 ! ynpp temperate conifer evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp05 = 162 ! ynpp temperate broadleaf cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp06 = 163 ! ynpp boreal conifer evergreen trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp07 = 164 ! ynpp boreal broadleaf cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp08 = 165 ! ynpp boreal conifer cold-deciduous trees
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp09 = 166 ! ynpp evergreen shrubs
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp10 = 167 ! ynpp cold-deciduous shrubs
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp11 = 168 ! ynpp warm (c4) grasses
  INTEGER, PUBLIC, PARAMETER :: nDiag_ynpp12 = 169 ! ynpp cool (c3) grasses

  INTEGER, PUBLIC, PARAMETER :: nDiag_cmontp = 170 !coldest monthly temperature                             (C)
  INTEGER, PUBLIC, PARAMETER :: nDiag_wmontp = 171 !warmest monthly temperature                             (C)
  INTEGER, PUBLIC, PARAMETER :: nDiag_atogpp = 172 !annual total gpp for ecosystem                          (kg-c/m**2/yr)
  INTEGER, PUBLIC, PARAMETER :: nDiag_toigpp = 173 !instantaneous gpp                                           (mol-CO2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_fxcsol = 174 !instantaneous fine co2 flux from soil                  (mol-CO2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_mcsoil = 175 !instantaneous microbial co2 flux from soil            (mol-CO2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cagcub = 176 !canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cagcuc = 177 !canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cagcls = 178 !canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cagcl4 = 179 !canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cagcl3 = 180 !canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cancub = 181 !canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cancuc = 182 !canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cancls = 183 !canopy average net photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cancl4 = 184 !canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cancl3 = 185 !canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cicoub = 186 !intercellular co2 concentration - broadleaf                  (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cicouc = 187 !intercellular co2 concentration - conifer                   (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cscoub = 188 !leaf boundary layer co2 concentration - broadleaf     (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_gscoub = 189 !upper canopy stomatal conductance - broadleaf          (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cscouc = 190 !leaf boundary layer co2 concentration - conifer          (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_gscouc = 191 !upper canopy stomatal conductance - conifer                  (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cicols = 192 !intercellular co2 concentration - shrubs                    (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cicol3 = 193 !intercellular co2 concentration - c3 plants                  (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cicol4 = 194 !intercellular co2 concentration - c4 plants                  (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cscols = 195 !leaf boundary layer co2 concentration - shrubs          (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_gscols = 196 !lower canopy stomatal conductance - shrubs                  (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cscol3 = 197 !leaf boundary layer co2 concentration - c3 plants     (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_gscol3 = 198 !lower canopy stomatal conductance - c3 grasses          (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_cscol4 = 199 !leaf boundary layer co2 concentration - c4 plants     (mol_co2/mol_air)
  INTEGER, PUBLIC, PARAMETER :: nDiag_gscol4 = 200 !lower canopy stomatal conductance - c4 grasses          (mol_co2 m-2 s-1)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tcthis = 201 !coldest monthly temperature of current year              (C)
  INTEGER, PUBLIC, PARAMETER :: nDiag_twthis = 202 !warmest monthly temperature of current year              (C)
 
  INTEGER, PUBLIC, PARAMETER :: nDiag_ObuLen = 203! Obukhov length                             (m)
  INTEGER, PUBLIC, PARAMETER :: nDiag_InPhiM = 204! inverse phi function for momentum          ()
  INTEGER, PUBLIC, PARAMETER :: nDiag_InPhiH = 205! inverse phi function for heat              ()
  INTEGER, PUBLIC, PARAMETER :: nDiag_Bouyac = 206! buoyancy scale (m/s**2)            
  INTEGER, PUBLIC, PARAMETER :: nDiag_qlicld = 207! liquid water content in cloud 
 
  INTEGER, PUBLIC, PARAMETER :: nDiag_shfcan = 208! sensible heat flux from canopy     (W m-2)
  INTEGER, PUBLIC, PARAMETER :: nDiag_shfgnd = 209! sensible heat flux from ground     (W m-2)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tracan = 210! transpiration from canopy          (W m-2)
  INTEGER, PUBLIC, PARAMETER :: nDiag_tragcv = 211! transpiration from ground cover    (W m-2)
  INTEGER, PUBLIC, PARAMETER :: nDiag_inlocp = 212! interception loss from canopy      (W m-2)
  INTEGER, PUBLIC, PARAMETER :: nDiag_inlogc = 213! interception loss from ground cover(W m-2)

  INTEGER, PUBLIC, PARAMETER :: nDiag_trcliq = 214 ! Liquid Mixing Ratio kg/kg
  INTEGER, PUBLIC, PARAMETER :: nDiag_tkemyj = 215 ! Turbulent Kinetic Energy
  INTEGER, PUBLIC, PARAMETER :: nDiag_trcice = 216 ! Ice Mixing Ratio kg/kg

  INTEGER, PUBLIC, PARAMETER :: nDiag_cape2d = 217 ! CONVECTIVE AVAIL. POT.ENERGY M2/S2	
  INTEGER, PUBLIC, PARAMETER :: nDiag_cine2d = 218 ! CONVECTIVE INHIB. ENERGY M2/S2
  INTEGER, PUBLIC, PARAMETER :: nDiag_sweath = 219 ! SEVERE WEATHER THREAT

  LOGICAL(KIND=i8), PUBLIC :: StartStorDiag = .FALSE. ! Start Storage Diagnostic

!--------------------------------------------------------------------------
  INTEGER                        :: nMax
  INTEGER                        :: mMax
  INTEGER                        :: mnMax
  INTEGER          , ALLOCATABLE :: mnMap(:,:)
  INTEGER          , ALLOCATABLE :: ibMaxPerJB(:)
  INTEGER                        :: iMaxNew
  INTEGER                        :: jMaxNew
  INTEGER                        :: kMaxNew
  INTEGER                        :: ibMax
  INTEGER                        :: jbMax
  INTEGER                        :: mxgaus
  INTEGER                        :: mxspec
  LOGICAL(KIND=i8)                        :: pfbar
  LOGICAL(KIND=i8)                        :: dodyn
  INTEGER                        :: mxrq
  INTEGER                        :: mgaus
  INTEGER                        :: ngaus
  INTEGER                        :: mspec
  INTEGER                        :: nspec
  INTEGER                        :: ispec
  INTEGER                        :: igaus
  INTEGER                        :: mxavl
  INTEGER                        :: icf
  INTEGER                        :: nof
  INTEGER                        :: ihdim (2)
  INTEGER          , ALLOCATABLE :: nuc   (:)
  CHARACTER(LEN=4 )         , ALLOCATABLE :: alias (:)
  CHARACTER(LEN=40)         , ALLOCATABLE :: chrdo (:)
  INTEGER          , ALLOCATABLE :: krcf  (:)
  INTEGER          , ALLOCATABLE :: jrcf  (:)
  INTEGER          , ALLOCATABLE :: kravl (:)
  INTEGER          , ALLOCATABLE :: ixcf  (:)
  INTEGER          , ALLOCATABLE :: incf  (:)
  INTEGER          , ALLOCATABLE :: iclcd (:)
  INTEGER          , ALLOCATABLE :: ixavl (:)
  INTEGER          , ALLOCATABLE :: inavl (:)
  LOGICAL(KIND=i8)          , ALLOCATABLE :: dodia (:)
  INTEGER          , ALLOCATABLE :: lvavl (:)
  INTEGER          , ALLOCATABLE :: nuavl (:)
  INTEGER          , ALLOCATABLE :: itavl (:)
  INTEGER          , ALLOCATABLE :: iavrq (:)
  INTEGER          , ALLOCATABLE :: nurq  (:)
  INTEGER          , ALLOCATABLE :: lvcf  (:)
  INTEGER          , ALLOCATABLE :: nucf  (:)
  INTEGER          , ALLOCATABLE :: itcf  (:)
  LOGICAL(KIND=i8)          , ALLOCATABLE :: icfu  (:)
  INTEGER          , ALLOCATABLE :: ixucf (:)
  INTEGER          , ALLOCATABLE :: inucf (:)
  INTEGER          , ALLOCATABLE :: jrucf (:)
  CHARACTER(len=40)         , ALLOCATABLE :: avail (:)
  CHARACTER(len=40)         , ALLOCATABLE :: reqdg (:)
  CHARACTER(len=40)         , ALLOCATABLE :: combf (:)
  REAL(KIND=r8)             , ALLOCATABLE :: gaus  (:,:,:)
  REAL(KIND=r8)             , ALLOCATABLE :: gaus_in(:,:,:)
  REAL(KIND=r8)             , ALLOCATABLE :: spec  (:,:)
  REAL(KIND=r8)             , ALLOCATABLE :: spec1d(:,:,:)
  INTEGER          , ALLOCATABLE :: lspec (:)
  INTEGER          , ALLOCATABLE :: lgaus (:)
  CHARACTER(LEN=16)         , ALLOCATABLE :: aunits(:) 
  REAL(KIND=r8)             , ALLOCATABLE :: dcol  (:)
  REAL(KIND=r8)             , ALLOCATABLE :: scol  (:)
  INTEGER          , ALLOCATABLE :: lvrq  (:)
  REAL(KIND=r8)             , ALLOCATABLE :: CountTOTAL(:,:)
  REAL(KIND=r8)             , ALLOCATABLE :: CountGaus  (:,:,:)
   
  INTERFACE updia
     MODULE PROCEDURE  updiaA,updiaB,StoreMaskedDiag1D,StoreMaskedDiag2D
  END INTERFACE

CONTAINS


  SUBROUTINE InitDiagnostics(mgaus_in  , ngaus_in , mspec_in, nspec_in, dodyn_in , colrad    , &
                             mMax_in   , nMax_in  , mnMax_in, mnMap_in,iMaxNew_in,jMaxNew_in , &
                             kMaxNew_in,ibMax_in  , jbMax_in,ibMaxPerJB_in, grid,fNameDTable ,&
                             fNameRDT)
    !
    !
    ! indiag :initialize diagnostic database; extended diagnostics version 1;
    !         this routine reads in the available diagnostics table and
    !         the desired diagnostic table; the two are compared and
    !         diagnostic tables are determined; these tables along with the
    !         standard prognostic output table are used to form the output
    !         directory; the actual accumulators are set in subroutine setdia;
    !         available diagnostics table should only be changed by appending;
    !         positions in table are permanently set when determined.
    !
    ! development notes
    !  
    !        version 1 of the extended diagnostics system removes data
    !        management of the diagnostic accumulators from the model code
    !        and allows for user selectable diagnostics.  this version will
    !        maintain the diagnostic accumulators in memory and permit the
    !        user to select individual diagnostics or to combine several into
    !        one diagnostic.  only available diagnostics can be selected or
    !        combined.
    !        later versions will allow for the use of solid state disk (ssd),
    !        regular disk, or other media to retain the accumulators.  this
    !        will allow for a large set of diagnostics with a reduced use of
    !        memory.  other changes will include use of dynamic memory and a
    !        user friendly interface.
    !
    
    INTEGER, INTENT(IN)    :: mgaus_in
    INTEGER, INTENT(IN)    :: ngaus_in
    INTEGER, INTENT(IN)    :: mspec_in
    INTEGER, INTENT(IN)    :: nspec_in
    LOGICAL(KIND=i8), INTENT(IN)    :: dodyn_in
    INTEGER, INTENT(IN)    :: mMax_in  
    INTEGER, INTENT(IN)    :: nMax_in  
    INTEGER, INTENT(IN)    :: mnMax_in 
    INTEGER, INTENT(IN)    :: mnMap_in(:,:)
    INTEGER, INTENT(IN)    :: iMaxNew_in
    INTEGER, INTENT(IN)    :: jMaxNew_in
    INTEGER, INTENT(IN)    :: kMaxNew_in
    INTEGER, INTENT(IN)    :: ibMax_in
    INTEGER, INTENT(IN)    :: jbMax_in
    INTEGER, INTENT(IN)    :: ibMaxPerJB_in(:)
    REAL(KIND=r8)   , INTENT(IN)    :: colrad(jMaxNew_in)
    CHARACTER(len=*), INTENT(IN)    :: grid    
    CHARACTER(len=*), INTENT(IN)    :: fNameDTable  
    CHARACTER(len=*), INTENT(IN)    :: fNameRDT

    INTEGER, ALLOCATABLE   :: jpavl(:)
    INTEGER, ALLOCATABLE   :: irqav(:)
    INTEGER, ALLOCATABLE   :: irqcf(:)
    LOGICAL(KIND=i8), ALLOCATABLE   :: irqu(:)
    INTEGER, ALLOCATABLE   :: kfrq(:)
    INTEGER, ALLOCATABLE   :: jpcf(:)
    CHARACTER(len=40)               :: ocf
    CHARACTER(len= 8)               :: typcd(2)
    CHARACTER(len=38)               :: poscd(0:3)
    INTEGER                :: diag
    INTEGER                :: ele
    INTEGER                :: l
    INTEGER                :: m
    INTEGER                :: n
    INTEGER                :: j
    INTEGER                :: iac
    INTEGER                :: nn
    INTEGER                :: ix
    INTEGER                :: ja
    INTEGER                :: ia
    INTEGER                :: k1
    INTEGER                :: k2
    INTEGER                :: k3
    INTEGER                :: i1
    INTEGER                :: i2
    INTEGER                :: mm
    INTEGER                :: kk
    INTEGER                :: jx
    INTEGER                :: la
    INTEGER                :: ka
    INTEGER                :: kka
    INTEGER                :: irix
    INTEGER                :: kx
    INTEGER                :: nofp
    INTEGER                :: mx
    INTEGER                :: in 
    REAL(KIND=r8)                   :: pie
    REAL(KIND=r8)                   :: colb(jMaxNew_in)
    
    mMax    = mMax_in
    nMax    = nMax_in
    mnMax   = mnMax_in   
    iMaxNew = iMaxNew_in
    jMaxNew = jMaxNew_in
    kMaxNew = kMaxNew_in
    ibMax   = ibMax_in
    jbMax   = jbMax_in
    sMapSize = nMax-1
    ALLOCATE(SMap(sMapSize))
    ALLOCATE (mnMap(mMax_in,nMax_in))
    ALLOCATE (ibMaxPerJB(jbMax))
    ibMaxPerJB=ibMaxPerJB_in
    mnMap = mnMap_in
    
    DO n = 1, nMax-1
       SMap(n) = 2*mnMap(1,n+1)-1
    END DO

    aMapSize = mnMax - nMax
    ALLOCATE(AMap(aMapSize))
    l = 0
    DO diag = 1, nMax 
       DO ele = 2, mMax-diag+1
          m = ele
          n = m+diag-1
          l = l+1
          AMap(l) = 2*mnMap(m,n)
       END DO
    END DO

    !
    !     avail = name of available diagnostic
    !     lvavl = levels in available diagnostic (1 or kmax)
    !     nuavl = unit code of available diagnostic
    !     itavl = type of available diagnostic (1 gaussian, 2 spectral)
    !     jpavl = position in code of available diagnostic (1 gloop/gfidi,
    !             2 gwater, 3 both, 0 neither)

    mgaus = mgaus_in
    ngaus = ngaus_in
    mspec = mspec_in
    nspec = nspec_in
    dodyn = dodyn_in




    ihdim(1)=iMaxNew*jMaxNew
    ihdim(2)=2*mnMax
    ihdim(2)=iMaxNew*jMaxNew
    ALLOCATE(dcol(jMaxNew))
    ALLOCATE(scol(jMaxNew))
    pie = 4.0_r8*ATAN(1.0_r8)
    !     
    !     define latitude grid for integration
    !     
    IF(grid == '2D  ') THEN
    colb(1) = 0.0_r8
    DO j = 2, jMaxNew
       colb(j) = 0.5_r8*(colrad(j)+colrad(j-1))
    END DO
    DO j = 1, jMaxNew-1
       dcol(j) = colb(j+1)-colb(j)
    END DO
    dcol(jMaxNew) = pie-colb(jMaxNew)
    DO j = 1, jMaxNew
       scol(j) = SIN(colrad(j))
    END DO
    ELSE
    colb(1) = 0.0_r8
    DO j = 1, jMaxNew
       colb(j) = 0.5_r8*(colrad(j))
    END DO
    DO j = 1, jMaxNew
       dcol(j) = colb(j)
    END DO
    dcol(jMaxNew) = pie-colb(jMaxNew)
    DO j = 1, jMaxNew
       scol(j) = SIN(colrad(j))
    END DO    
    END IF
    ALLOCATE(aunits(-1:numx))
    ALLOCATE(dodia(ndavl))
    ALLOCATE(lvavl(ndavl))
    ALLOCATE(nuavl(ndavl))
    ALLOCATE(itavl(ndavl))
    ALLOCATE(iavrq(ndavl))
    ALLOCATE(ixavl(ndavl))
    ALLOCATE(inavl(ndavl))
    ALLOCATE(nurq (ndrq ))
    ALLOCATE(iclcd(ndrq ))
    ALLOCATE(lvcf (ncdg ))
    ALLOCATE(nucf (ncdg ))
    ALLOCATE(ixcf (ncdg ))
    ALLOCATE(incf (ncdg ))
    ALLOCATE(itcf (ncdg ))
    ALLOCATE(icfu (ncdg ))
    ALLOCATE(ixucf(ncdg ))
    ALLOCATE(inucf(ncdg ))
    ALLOCATE(kravl(jxavl))
    ALLOCATE(krcf (jxcdg))
    ALLOCATE(jrcf (jxcdg))
    ALLOCATE(jrucf(jxcdg))  
    ALLOCATE(avail(ndavl))
    ALLOCATE(reqdg(ndrq) )
    ALLOCATE(combf(ncdg) )
    ALLOCATE(lspec(-ncdg:ndavl))
    ALLOCATE(lgaus(-ncdg:ndavl))


    ALLOCATE(jpavl(ndavl))
    ALLOCATE(irqav(ndrq ))
    ALLOCATE(irqcf(ndrq ))
    ALLOCATE(irqu (ndrq ))
    ALLOCATE(lvrq (ndrq ))
    ALLOCATE(kfrq (ncdg )) 
    ALLOCATE(jpcf (ncdg ))

    ALLOCATE(nuc  (ndrq))
    ALLOCATE(alias(ndrq))
    ALLOCATE(chrdo(ndrq))

    combf=' '
    jpcf=0


  ! field name     !     avail = name of available diagnostic
 
  avail="                                        "
    avail(1:39)=(/  &
         'TIME MEAN SURFACE PRESSURE              ', &
         'TIME MEAN DIVERGENCE                    ', &
         'TIME MEAN VORTICITY                     ', &
         'TIME MEAN SPECIFIC HUMIDITY             ', &
         'TIME MEAN VIRTUAL TEMPERATURE           ', &
         'TIME MEAN SURFACE TEMPERATURE           ', &
         'TIME MEAN OMEGA                         ', &
         'TIME MEAN SIGMADOT                      ', &
         'TOTAL PRECIPITATION                     ', &
         'CONVECTIVE PRECIPITATION                ', &
         'LARGE SCALE PRECIPITATION               ', &
         'SNOWFALL                                ', &
         'RUNOFF                                  ', &
         'PRECIPITABLE WATER                      ', &
         'INTERCEPTION LOSS                       ', &
         'SENSIBLE HEAT FLUX FROM SURFACE         ', &
         'LATENT HEAT FLUX FROM SURFACE           ', &
         'SURFACE ZONAL WIND STRESS               ', &
         'SURFACE MERIDIONAL WIND STRESS          ', &
         'CLOUD COVER                             ', &
         'DOWNWARD LONG WAVE AT BOTTOM            ', &
         'UPWARD LONG WAVE AT BOTTOM              ', &
         'OUTGOING LONG WAVE AT TOP               ', &
         'INCIDENT SHORT WAVE FLUX                ', &
         'DOWNWARD SHORT WAVE AT GROUND           ', &
         'UPWARD SHORT WAVE AT GROUND             ', &
         'UPWARD SHORT WAVE AT TOP                ', &
         'SHORT WAVE ABSORBED BY EARTH/ATMOSPHERE ', &
         'SHORT WAVE ABSORBED AT GROUND           ', &
         'NET LONG WAVE AT BOTTOM                 ', &
         'LONG WAVE RADIATIVE HEATING             ', &
         'SHORT WAVE RADIATIVE HEATING            ', &
         'CONVECTIVE LATENT HEATING               ', &
         'CONVECTIVE MOISTURE SOURCE              ', &
         'LARGE SCALE LATENT HEATING              ', &
         'LARGE SCALE MOISTURE SOURCE             ', &
         'SHALLOW CONVECTIVE HEATING              ', &
         'SHALLOW CONV. MOISTURE SOURCE           ', &
         'VERTICAL DIFFUSION HEATING              '/)
    avail(40:78)=(/  &
         'VERTICAL DIFF. MOISTURE SOURCE          ', &
         'VERTICAL DIFFUSION DU/DT                ', &
         'VERTICAL DIFFUSION DV/DT                ', &
         'GRAVITY WAVE DRAG SFC ZONAL STRESS      ', &
         'GRAVITY WAVE DRAG SFC MERIDIONAL STRESS ', &
         'GRAVITY WAVE DRAG DU/DT                 ', &
         'GRAVITY WAVE DRAG DV/DT                 ', &
         'HORIZONTAL HEATING DIFFUSION            ', &
         'HORIZONTAL MOISTURE DIFFUSION           ', &
         'HORIZONTAL DIVERGENCE DIFFUSION         ', &
         'HORIZONTAL VORTICITY DIFFUSION          ', &
         'DIVERGENCE * SPECIFIC HUMIDITY          ', &
         'VERTICAL MOISTURE ADVECTION             ', &
         'HORIZ. MOISTURE FLUX CONV.              ', &
         'VERT. INTEGRATED MOISTURE FLUX CONV.    ', &
         'TIME MEAN LOG SURFACE PRESSURE          ', &
         'DOWNWARD LONG WAVE AT BOTTOM (CLEAR)    ', &
         'OUTGOING LONG WAVE AT TOP (CLEAR)       ', &
         'DOWNWARD SHORT WAVE AT GROUND (CLEAR)   ', &
         'UPWARD SHORT WAVE AT GROUND (CLEAR)     ', &
         'UPWARD SHORT WAVE AT TOP (CLEAR)        ', &
         'SHORT WV ABSRBD BY EARTH/ATMOS (CLEAR)  ', &
         'SHORT WAVE ABSORBED AT GROUND (CLEAR)   ', &
         'NET LONG WAVE AT BOTTOM (CLEAR)         ', &
         'TIME MEAN DEEP SOIL TEMPERATURE         ', &
         'GROUND/SURFACE COVER TEMPERATURE        ', &
         'CANOPY TEMPERATURE                      ', &
         'TEMPERATURE OF CANOPY AIR SPACE         ', &
         'VAPOR PRESSURE OF CANOPY AIR SPACE      ', &
         'BARE SOIL LATENT HEAT                   ', &
         'NEG. HUM. CORR. MOISTURE SOURCE         ', &
         'OZONE MIXING RATIO                      ', &
         'VERTICAL DIST TOTAL CLOUD COVER         ', &
         'INVERSION CLOUD                         ', &
         'SUPERSATURATION CLOUD                   ', &
         'CONVECTIVE CLOUD                        ', &
         'SHALLOW CONVECTIVE CLOUD                ', &
         'CLOUD LIQUID WATER PATH                 ', &
         'LONGWAVE CLOUD EMISSIVITY               '/)
    avail(79:133)=(/  &
         'SHORTWAVE CLOUD OPTICAL DEPTH           ', &
         'CANOPY AIR SPC TO REF. LVL RESISTANCE   ', &
         'CANOPY AIR SPC TO CANOPY RESISTANCE     ', &
         'CANOPY AIR SPC TO GROUND RESISTANCE     ', &
         'CANOPY RESISTANCE                       ', &
         'GROUND COVER RESISTANCE                 ', &
         'BARE SOIL SURFACE RESISTANCE            ', &
         'HORIZONTAL MOMENTUM TRANSPORT           ', &
         'VERTICAL ZONAL MOMENTUM TRANSPORT       ', &
         'VERTICAL MERIDIONAL MOMENTUM TRANSPORT  ', &
         'MERIDIONAL SENSIBLE HEAT TRANSPORT      ', &
         'ZONAL SENSIBLE HEAT TRANSPORT           ', &
         'VERTICAL SENSIBLE HEAT TRANSPORT        ', &
         'MERIDIONAL SPECIFIC HUMIDITY TRANSPORT  ', &
         'ZONAL SPECIFIC HUMIDITY TRANSPORT       ', &
         'VERTICAL SPECIFIC HUMIDITY TRANSPORT    ', &
         'LONG WAVE RADIATIVE HEATING (CLEAR)     ', & !hmjb
         'SHORT WAVE RADIATIVE HEATING (CLEAR)    ', & !hmjb
         'TIME MEAN TEMP AT 2-M FROM SFC          ', &
         'TIME MEAN SPEC HUMIDITY AT 2-M FROM SFC ', &
         'SPEED WIND AT 2-M FROM SURFACE          ', &
         'TIME MEAN TEMP AT 10-M FROM SFC         ', &
         'TIME MEAN SPEC HUMIDITY AT 10-M FROM SFC', &
         'SPEED WIND AT 10-M FROM SURFACE         ', &
         'VERTICALLY INTEGRATED OZONE CONTENT     ', &!hmjb
         'DEW POINT TEMPERATURE                   ', &
         'TIME MEAN AT 10 METRE U-WIND COMPONENT  ', &
         'TIME MEAN AT 10 METRE V-WIND COMPONENT  ', &
         'TOT BIOMASS IN THE UPPER CANOPY         ', &
         'TOT BIOMASS IN THE LOWER CANOPY         ', &
         'TOT LEAF AREA INDEX FOR THE UPPER CANOPY', &
         'TOT LEAF AREA INDEX FOR THE LOWER CANOPY', &
         'TOT STORAGE OF N IN SOIL PROFILE        ', &
         'TOT. OF WATER STORAGE SNOW SOIL VEG     ', &
         'LITTER DECOMPOSITION FACTOR             ', &
         'SOIL ORGANIC MATTER DECOMPOSITION FACTOR', &
         'FRAC OVERALL AREA COVER BY UPPER CANOPY ', &
         'FRAC SNOWFREE AREA COVER BY LOWER CANOPY', &
         'FRACTIONAL SNOW COVER                   ', &
         'INSTANTANEOUS NPP                       ', &
         'INS. NET ECOSY. EXCHANGE CO2 P TIMESTEP ', &
         'ANNUAL TOT GROW DEGREE DAYS > 0C        ', &
         'ANNUAL TOT GROW DEGREE DAYS > 5C        ', &
         'MONTH AVE 2-M SURFACE-AIR TEMPERATURE   ', &
         'MONTHLY TOTAL NPP FOR ECOSYSTEM         ', &
         'MONTH TOT NET ECOSYSTEM EXCHANGE OF CO2 ', &
         'ANNUAL TOTAL NPP FOR ECOSYSTEM          ', &
         'ANNUAL TOTAL NEE FOR ECOSYSTEM          ', &
         'UPPER CANOPY SINGLE-SIDED LAI           ', &
         'LOWER CANOPY SINGLE-SIDED LAI           ', &
         'SURFACE FRICTION VELOCITY               ', &
         'PLANETARY BOUNDARY LAYER HEIGHT         ', &
         'DIFFUSION COEFFICIENT FOR HEAT          ', &
         'DIFFUSION COEFFICIENT FOR MOMENTUM      ', &
         'BULK RICHARDSON NO. REF LEVEL           '/)
   avail(134:169)=(/  &
         'PFT TROPICAL BROADLEAF EVERGREEN TREES  ', &
         'PFT TR. BROADLEAF DROUGHT-DECIDUOUS T.  ', &
         'PFT WARM-TEMPERATE BROADLEAF EVERGREEN  ', &
         'PFT TEMPERATE CONIFER EVERGREEN TREES   ', & 
         'PFT TEMPERATE BROADLEAF COLD-DECIDUOUS  ', &
         'PFT BOREAL CONIFER EVERGREEN TREES      ', &
         'PFT BOREAL BROADLEAF COLD-DECIDUOUS     ', &
         'PFT BOREAL CONIFER COLD-DECIDUOUS       ', &
         'PFT EVERGREEN SHRUBS                    ', &
         'PFT COLD-DECIDUOUS SHRUBS               ', &
         'PFT WARM (C4) GRASSES                   ', &
         'PFT COOL (C3) GRASSES                   ', &
         'CBIOL TROPICAL BROADLEAF EVERGREEN      ', &
         'CBIOL TR. BROADLEAF DROUGHT-DECIDUOUS   ', &
         'CBIOL WARM-TEMPERATE BROADLEAF EV.GREEN ', &
         'CBIOL TEMPERATE CONIFER EVERGREEN TREES ', &
         'CBIOL TEMPERATE BRD.LEAF COLD-DECIDUOUS ', &
         'CBIOL BOREAL CONIFER EVERGREEN TREES    ', &
         'CBIOL BOREAL BROADLEAF COLD-DECIDUOUS   ', &
         'CBIOL BOREAL CONIFER COLD-DECIDUOUS     ', &
         'CBIOL EVERGREEN SHRUBS                  ', &
         'CBIOL COLD-DECIDUOUS SHRUBS             ', &
         'CBIOL WARM (C4) GRASSES                 ', &
         'CBIOL COOL (C3) GRASSES                 ', &
         'YNPP TROPICAL BROADLEAF EVERGREEN TREES ', &
         'YNPP TR BROADLEAF DROUGHT-DECIDUOUS T.  ', &
         'YNPP WARM-TEMPERATE BROADLEAF EVERGREEN ', &
         'YNPP TEMPERATE CONIFER EVERGREEN TREES  ', &
         'YNPP TEMPERATE BROADLEAF COLD-DECIDUOUS ', &
         'YNPP BOREAL CONIFER EVERGREEN TREES     ', &
         'YNPP BOREAL BROADLEAF COLD-DECIDUOUS    ', &
         'YNPP BOREAL CONIFER COLD-DECIDUOUS      ', &
         'YNPP EVERGREEN SHRUBS                   ', &
         'YNPP COLD-DECIDUOUS SHRUBS              ', &
         'YNPP WARM (C4) GRASSES                  ', &
         'YNPP COOL (C3) GRASSES                  '/)
   avail(170:ndavl)=(/  &
         'COLDEST MONTHLY TEMPERATURE             ', &
         'WARMEST MONTHLY TEMPERATURE             ', &
         'ANNUAL TOTAL GPP FOR ECOSYSTEM          ', &
         'INSTANTANEOUS GPP                       ', &
         'INSTANTANEOUS FINE CO2 FLUX FROM SOIL   ', &
         'INSTAN.  MICROBIAL CO2 FLUX FROM SOIL   ', &
         'GROSS PHOTOSYNTHESIS RATE - BROADLEAF   ', &
         'GROSS PHOTOSYNTHESIS RATE - CONIFER     ', &
         'GROSS PHOTOSYNTHESIS RATE - SHRUBS      ', &
         'GROSS PHOTOSYNTHESIS RATE - C4 GRASSES  ', &
         'GROSS PHOTOSYNTHESIS RATE - C3 GRASSES  ', &
         'NET PHOTOSYNTHESIS RATE - BROADLEAF     ', &
         'NET PHOTOSYNTHESIS RATE - CONIFER       ', &
         'NET PHOTOSYNTHESIS RATE - SHRUBS        ', &
         'NET PHOTOSYNTHESIS RATE - C4 GRASSES    ', &
         'NET PHOTOSYNTHESIS RATE - C3 GRASSES    ', &
         'INTERCELLULAR CO2CONCENTRATION BROADLEAF', &
         'INTERCELLULAR CO2CONCENTRATION CONIFER  ', &
         'LEAF BL CO2 CONCENTRATION - BROADLEAF   ', &
         'UPPER CAS STOMATAL CONDUCTANCE BROADLEAF', &
         'LEAF BL CO2 CONCENTRATION - CONIFER     ', &
         'UPPER CAS STOMATAL CONDUCTANCE - CONIFER', &
         'INTERCELLULAR CO2CONCENTRATION SHRUBS   ', &
         'INTERCELLULAR CO2CONCENTRATION C3 PLANTS', &
         'INTERCELLULAR CO2CONCENTRATION C4 PLANTS', &
         'LEAF BL CO2 CONCENTRATION - SHRUBS      ', &
         'LOWER CAS STOMATAL CONDUCTANCE - SHRUBS ', &
         'LEAF BL CO2 CONCENTRATION - C3 PLANTS   ', &
         'LOWER CAS STOMATAL CONDUCTANCE C3GRASSES', &
         'LEAF BL CO2 CONCENTRATION - C4 PLANTS   ', &
         'LOWER CAS STOMATAL CONDUCTANCE C4GRASSES', &
         'COLDEST MON TEMPERATURE OF CURRENT YEAR ', &
         'WARMEST MON TEMPERATURE OF CURRENT YEAR ', &
         'OBUKHOV LENGTH                          ', &
         'INVERSE PHI FUNCTION FOR MOMENTUM       ', &
         'INVERSE PHI FUNCTION FOR HEAT           ', &
         'BUOYANCY SCALE                          ', &
         'LIQUID WATER CONTENT IN CLOUD           ', &
         'SENSIBLE HEAT FLUX FROM CANOPY          ', &
         'SENSIBLE HEAT FLUX FROM GROUND          ', &
         'TRANSPIRATION FROM CANOPY               ', &
         'TRANSPIRATION FROM GROUND COVER         ', &
         'INTERCEPTION LOSS FROM CANOPY           ', &
         'INTERCEPTION LOSS FROM GROUND COVER     ', &
         'LIQUID MIXING RATIO KG/KG               ', &
         'TURBULENT KINETIC ENERGY                ', &
         'ICE MIXING RATIO KG/KG                  ', &
         'CONVECTIVE AVAIL. POT.ENERGY            ', &
         'CONVECTIVE INHIB. ENERGY                ', &
         'SEVERE WEATHER THREAT                   '/)

    !     lvavl = levels in available diagnostic (1 or 2=kmax)

    lvavl(1:ndavl)=(/ &
         1,    2,    2,    2,    2,    1,    2,    2,    1,    1, &!   1  -  10
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &!  11  -  20
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &!  21  -  30
         2,    2,    2,    2,    2,    2,    2,    2,    2,    2, &!  31  -  40
         2,    2,    1,    1,    2,    2,    2,    2,    2,    2, &!  41  -  50
         2,    2,    2,    1,    1,    1,    1,    1,    1,    1, &!  51  -  60
         1,    1,    1,    1,    1,    1,    1,    1,    1,    2, &!  61  -  70
         2,    2,    2,    2,    2,    2,    2,    2,    2,    1, &!  71  -  80
         1,    1,    1,    1,    1,    2,    2,    2,    2,    2, &!  81  -  90
         2,    2,    2,    2,    2,    2,    1,    1,    1,    1, &!  91  - 100
         1,    1,    1,    2,    1,    1,    1,    1,    1,    1, &! 101  - 110
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 111  - 120
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 121  - 130
         2,    2,    2,    1,    1,    1,    1,    1,    1,    1, &! 131  - 140
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 141  - 150
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 151  - 160
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 161  - 170
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 171  - 180
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 181  - 190
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &! 191  - 200
         1,    1,    1,    1,    1,    1,    2,    1,    1,    1, &! 201  - 210
         1,    1,    1,    2,    2,    2,    1,    1,    1/)

    !     nuavl = unit code of available diagnostic

    nuavl(1:ndavl)=(/ &
         132,  50,  50,   0,  40,  40, 153,  50, 120, 120, &!   1  -  10
         120, 120, 120, 110, 170, 170, 170, 130, 130,   0, &!  11  -  20
         170, 170, 170, 170, 170, 170, 170, 170, 170, 170, &!  21  -  30
         70 ,  70,  70,  50,  70,  50,  70,  50,  70,  50, &!  31  -  40
         100, 100, 130, 130, 100, 100,  70,  50,  80,  80, &!  41  -  50
         50 ,  50,  50, 120, 142, 170, 170, 170, 170, 170, &!  51  -  60
         170, 170, 170,  40,  40,  40,  40, 131, 170,  50, &!  61  -  70
         0  ,   0,   0,   0,   0,   0,   0,   0,   0, 190, &!  71  -  80
         190, 190, 190, 190, 190, 180, 252, 252, 230, 230, &!  81  -  90
         242,  60,  60, 153,  70,  70,  40,   0,  60,  40, &!  91  - 100
         0  ,  60,   0,  40,  60,  60, 260, 260,   0,   0, &! 101  - 110
         261, 262,   0,   0,   1,   1,   1, 280,   0,  41, &! 111  - 120
         41 ,  40, 271, 271, 270, 270,   0,   0,  60,  10, &! 121  - 130
         90 ,  90,   0,   0,   0,   0,   0,   0,   0,   0, &! 131  - 140
         0,     0,   0,   0,   0, 260, 260, 260, 260, 260, &! 141  - 150
         260, 260, 260, 260, 260, 260, 260, 270, 270, 270, &! 151  - 160
         270, 270, 270, 270, 270, 270, 270, 270, 270,  41, &! 161  - 170
         41 , 270, 280, 280, 280, 280, 280, 280, 280, 280, &! 171  - 180
         280, 280, 280, 280, 280, 281, 281, 281, 280, 281, &! 181  - 190
         280, 281, 281, 281, 281, 280, 281, 280, 281, 280, &! 191  - 200
         41 ,  41,  10,   0,   0, 100,   0, 170, 170, 170, &! 201  - 210
         170, 170, 170,   0,  60,   0, 180, 180,   0 /)

    !     itavl = type of available diagnostic (1 gaussian, 2 spectral)
        
    itavl(1:ndavl)=(/ &
         1,    2,    2,    2,    2,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    2,    2,    2,    2, &
         1,    1,    2,    2,    2,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,    1,    1,    1,    1,    1/)

    !     jpavl = position in code of available diagnostic (1 gloop/gfidi,
    !             2 gwater, 3 both, 0 neither)
    
    jpavl(1:ndavl)=(/ &
         0,    0,    0,    0,   0,    1,    1,    1,    2,    2, &
         2,    2,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    2,    2,   2,    2,    2,    2,    1,    1, &
         1,    1,    1,    1,   1,    1,    0,    0,    0,    0, &
         1,    1,    0,    0,   0,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    2, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1,    1, &
         1,    1,    1,    1,   1,    1,    1,    1,    1/)

    IF(ifprt(51).GE.1)WRITE(nfprt,30)
    mxavl=0
    DO m=1,ndavl 
       IF (avail(m)(1:5) .NE. "     ") THEN
         IF (lvavl(m) .NE. 1) lvavl(m)=kMaxNew
         IF(ifprt(51).GE.1)WRITE(nfprt,60)avail(m),lvavl(m),nuavl(m), &
                                          itavl(m),jpavl(m)
         dodia(m)=.FALSE.
         iavrq(m)=0
         mxavl=mxavl+1
       END IF
    END DO

    IF(mxavl.LE.0)THEN
       WRITE(nfprt,6600)
       WRITE(nferr,6600)
       STOP 6600
    END IF

    !     reqdg = name of requested diagnostic
  reqdg="                                        "
  reqdg(1:20)=(/ &
      "TOTAL PRECIPITATION                     ", &
      "CONVECTIVE PRECIPITATION                ", &
      "LARGE SCALE PRECIPITATION               ", &
      "SNOWFALL                                ", &
      "RUNOFF                                  ", &
      "INTERCEPTION LOSS                       ", &
      "SENSIBLE HEAT FLUX FROM SURFACE         ", &
      "LATENT HEAT FLUX FROM SURFACE           ", & 
      "SURFACE ZONAL WIND STRESS               ", &
      "SURFACE MERIDIONAL WIND STRESS          ", &
      "CLOUD COVER                             ", &
      "DOWNWARD LONG WAVE AT BOTTOM            ", &
      "UPWARD LONG WAVE AT BOTTOM              ", &
      "OUTGOING LONG WAVE AT TOP               ", &
      "DOWNWARD SHORT WAVE AT GROUND           ", &
      "UPWARD SHORT WAVE AT GROUND             ", &
      "UPWARD SHORT WAVE AT TOP                ", &
      "SHORT WAVE ABSORBED AT GROUND           ", &
      "NET LONG WAVE AT BOTTOM                 ", &
      "GROUND/SURFACE COVER TEMPERATURE        "/)

    !     lvrq  = levels in requested diagnostic (1 or kmax)

  lvrq(1:20)=(/ &
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1, &
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1 /)
    
    !     nurq  = unit code of requested diagnostic

  nurq(1:20)=(/ &
  121,  121,  121,  121,  121,  170,  170,  170,  130,  130, &
    0,  170,  170,  170,  170,  170,  170,  170,  170,   40/)
    
    !     iclcd = requested diagnostic calculation code (0 direct
    !             calculation, > 0 add to requested field number iclcd,
    !             < 0 subtract from requested field number -iclcd )

  iclcd(1:20)=(/ &
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, &
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0/)

    IF(ifprt(51).GE.1)WRITE(nfprt,130)
    OPEN(49, file=TRIM(fNameDTable),ACTION="read")  
    mxrq=1
    DO n=1,ndrq
       READ(49,'(A40,3I5)',END=225)reqdg(n),lvrq(n),nurq(n),iclcd(n)
       IF (reqdg(n)(1:5) .NE. "     ")THEN
         IF (lvrq(n) .NE. 1) lvrq(n)=kMaxNew
         IF(ifprt(51).GE.1)WRITE(nfprt,160)reqdg(n),lvrq(n),nurq(n), &
                                           iclcd(n)
         irqcf(n)=0
         irqav(n)=0
         irqu(n)=.FALSE.
         mxrq=mxrq+1
       END IF 
    END DO
225 mxrq=n-1
    CLOSE(49,STATUS='KEEP') 
    OPEN(132, file=TRIM(fNameRDT),ACTION="read")   
    DO n = 1, mxrq
       READ(132,'(a40,i5,1x,a4)',END=40) chrdo(n),nuc(n),alias(n)
    END DO
40  n=0   
    CLOSE(132,STATUS='KEEP')   
    IF(mxrq.LE.0)THEN
       WRITE(nfprt,7100)
       WRITE(nferr,7100)
       STOP 7100
    END IF

    typcd(1)='GAUSSIAN'
    typcd(2)='SPECTRAL'
    poscd(0)='NOT COMPUTED IN EITHER GFIDI OR GWATER'
    poscd(1)='COMPUTED ONLY IN GFIDI                '
    poscd(2)='COMPUTED ONLY IN GWATER               '
    poscd(3)='COMPUTED IN BOTH GFIDI AND GWATER     '
    !     search for combined field components.  save as combined fields
    !     those fields which have at least one component.  mark desired
    !     field refered by component as a valid combined field
    !     (irqcf=-icf).  indicate combined field for component (irqcf=icf).
    nof=999999
    icf=0

    DO n=1,mxrq

       IF(nof.EQ.999999)THEN
          IF(iclcd(n).NE.0)nof=n-1
       END IF

       IF(nof.NE.999999)THEN

          IF(iclcd(n).EQ.0)THEN
             WRITE(nfprt,2100)n,nof
             WRITE(nferr,2100)n,nof
             STOP 2100
          END IF

          iac=abs(iclcd(n))

          IF(iac.GE.n)THEN
             WRITE(nfprt,2600)iac,n
             WRITE(nferr,2600)iac,n
             STOP 2600
          END IF

          IF(iclcd(iac).NE.0)THEN
             WRITE(nfprt,3100)iac,n,iclcd(iac)
             WRITE(nferr,3100)iac,n,iclcd(iac)
             STOP 3100
          END IF

          IF(reqdg(iac).EQ.reqdg(n))THEN
             WRITE(nfprt,3600)iac,n,reqdg(n)
             WRITE(nferr,3600)iac,n,reqdg(n)
             STOP 3600
          END IF

          irqu(n)=.TRUE.

          IF(icf.NE.0)THEN
             DO ix=1,icf
                IF(reqdg(iac).EQ.combf(ix))go to 270
             END DO
          END IF

          icf=icf+1
          combf(icf)=reqdg(iac)
          lvcf(icf)=lvrq(iac)
          nucf(icf)=nurq(iac)
          kfrq(icf)=iac
          irqcf(iac)=-icf
          irqu(iac)=.TRUE.
270       irqcf(n)=icf
       END IF

    END DO

    IF(nof.EQ.999999) nof=mxrq
    !     find available diagnostics corresponding to desired diagnostic.
    !     first, find directly available desired diagnostic
    DO nn=1,nof

       IF(irqcf(nn).NE.0)CYCLE

       DO m=1,mxavl
          IF(reqdg(nn).EQ.avail(m))go to 360
       END DO

       WRITE(nfprt,4100)nn,reqdg(nn)
       WRITE(nferr,4100)nn,reqdg(nn)
       STOP 4100

360    irqu(nn)=.TRUE.
       dodia(m)=.TRUE.
       irqav(nn)=m
       iavrq(m)=nn

    END DO
    !     second, find available diagnostic components for combined fields
    !     or find other combined fields used as components.  save
    !     component index (+ for a.d., - for c.f.)
    IF(nof.LT.mxrq)THEN
       nofp=nof+1
       DO nn=nofp,mxrq

          DO m=1,mxavl
             IF(reqdg(nn).EQ.avail(m))go to 480
          END DO

          IF(icf.NE.0)THEN
             DO ix=1,icf
                IF(reqdg(nn).EQ.combf(ix))go to 490
             END DO
          END IF

          WRITE(nfprt,4600)nn,reqdg(nn)
          WRITE(nferr,4600)nn,reqdg(nn)
          STOP 4600

480       dodia(m)=.TRUE.
          irqav(nn)=m
          go to 495

490       irqav(nn)=-ix
495       irqu(nn)=.TRUE.

       END DO
    END IF

    !     check to make sure all desired diagnostics are used

    DO n=1,mxrq
       IF(.NOT.irqu(n))THEN

          WRITE(nfprt,5100)n,reqdg(n)
          WRITE(nferr,5100)n,reqdg(n)
          STOP 5100

       END IF
    END DO

    !     find all components for each combined field

    IF(icf.NE.0)THEN
       ja=1
       DO ix=1,icf
          ixcf(ix)=ja
          itcf(ix)=0
          ia=0
          k1=0
          k2=0
          k3=0
          i1=0
          i2=0
          DO n=1,mxrq

             IF(irqcf(n).NE.ix)CYCLE
             ia=ia+1
             jrcf(ja)=n

             IF(irqav(n).GT.0)THEN

                !     case for available diagnostic component

                mm=irqav(n)

                IF(nuavl(mm)/10.NE.nucf(ix)/10)THEN
                   WRITE(nfprt,9100)nuavl(mm),mm,nucf(ix),ix
                   WRITE(nferr,9100)nuavl(mm),mm,nucf(ix),ix
                   STOP 9100
                END IF

                krcf(ja)=mm
                IF(jpavl(mm).EQ.1)k1=k1+1
                IF(jpavl(mm).EQ.2)k2=k2+1
                IF(jpavl(mm).EQ.3)k3=k3+1
                IF(itavl(mm).EQ.1)i1=i1+1
                IF(itavl(mm).EQ.2)i2=i2+1

             ELSE IF(irqav(n).LT.0)THEN

                !     case for combined field component
                kk=-irqav(n)

                IF(nucf(ix).NE.nucf(kk))THEN
                   WRITE(nfprt,9600)nucf(kk),kk,nucf(ix),ix
                   WRITE(nferr,9600)nucf(kk),kk,nucf(ix),ix
                   STOP 9600
                END IF

                krcf(ja)=irqav(n)
                IF(jpcf(kk).EQ.1)k1=k1+1
                IF(jpcf(kk).EQ.2)k2=k2+1
                IF(jpcf(kk).EQ.3)k3=k3+1
                IF(itcf(kk).EQ.1)i1=i1+1
                IF(itcf(kk).EQ.2)i2=i2+1

             END IF

             ja=ja+1

             IF(ja.GT.jxcdg)THEN
                WRITE(nfprt,5600)ix,n
                WRITE(nferr,5600)ix,n
                STOP 5600
             END IF

          END DO

          incf(ix)=ia

          IF(k3.NE.0)THEN
             jpcf(ix)=3
          ELSE IF(k2.NE.0.AND.k1.NE.0)THEN
             jpcf(ix)=3
          ELSE IF(k2.NE.0)THEN
             jpcf(ix)=2
          ELSE IF(k1.NE.0)THEN
             jpcf(ix)=1
          ELSE
             jpcf(ix)=0
          END IF

          IF(i1.NE.0.AND.i2.NE.0)THEN
             WRITE(nfprt,7600)ix,combf(ix),(krcf(mx),mx=ixcf(ix),ja-1)
             WRITE(nferr,7600)ix,combf(ix),(krcf(mx),mx=ixcf(ix),ja-1)
             STOP 7600
          END IF

          IF(i1.NE.0)itcf(ix)=1
          IF(i2.NE.0)itcf(ix)=2

       END DO
    END IF

    !     determine all available diagnostic uses

    ja=1
    DO m=1,mxavl

       IF(.NOT.dodia(m))CYCLE

       ixavl(m)=ja
       ia=0

       IF(iavrq(m).GT.0)THEN
          ia=ia+1
          kravl(ja)=iavrq(m)
          ja=ja+1

          IF(ja.GT.jxavl)THEN
             WRITE(nfprt,6100)m,ix
             WRITE(nferr,6100)m,ix
             STOP 6100
          END IF

       END IF

       IF(icf.NE.0)THEN
          DO ix=1,icf
             ka=ixcf(ix)
             DO in=1,incf(ix)
                IF(krcf(ka).EQ.m)THEN
                   ia=ia+1
                   kravl(ja)=-ix
                   ja=ja+1
                   IF(ja.GT.jxavl)THEN
                      WRITE(nfprt,6100)m,ix
                      WRITE(nferr,6100)m,ix
                      STOP 6100
                   END IF
                END IF
                ka=ka+1
             END DO
          END DO
       END IF

       inavl(m)=ia
    END DO
    !     find combined fields requiring given combined field component
    la=1
    DO ix=1,icf
       icfu(ix)=.FALSE.
       ixucf(ix)=0
       ia=0
       DO jx=1,icf
          ka=ixcf(jx)
          DO in=1,incf(jx)
             IF(krcf(ka).EQ.-ix)go to 840
             ka=ka+1
          END DO
          CYCLE
840       IF(jx.LE.ix)THEN
             WRITE(nfprt,8600)jx,ix,jrcf(jx),jrcf(ix)
             WRITE(nferr,8600)jx,ix,jrcf(jx),jrcf(ix)
             STOP 8600
          END IF
          IF(.NOT.icfu(ix))THEN
             icfu(ix)=.TRUE.
             ixucf(ix)=la
          END IF
          IF(iclcd(jrcf(jx)).GT.0)THEN
             jrucf(la)=jx
          ELSE
             jrucf(la)=-jx
          END IF
          ia=ia+1
          la=la+1
          IF(la.GT.jxcdg)THEN
             WRITE(nfprt,8100)jx,ix
             WRITE(nferr,8100)jx,ix
             STOP 8100
          END IF
       END DO
       inucf(ix)=ia
    END DO
    
    !     read in conversion tables
    
    aunits(:)=(/ &
         "Unknown         ","No Dim          ","%               ",&
         "Gm/Kg           ","Ppm             ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","M               ",&
         "Cm              ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Kg              ","Gm              ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Sec             ","Days            ",&
         "Yrs             ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","K               ",&
         "C               ","F               ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "1/Sec           ","1/Day           ","Gm/Kg/Day       ",&
         "10**-5 1/Sec    ","10**-6 1/Sec    ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","M/Sec           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","K/Sec           ",&
         "K/Day           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Sec**-2         ","1/Sec/Day       ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","M**2/Sec        ","10**6 M**2/Sec  ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","M Sec**-2       ",&
         "M/Sec/Day       ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Kg M**-2        ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Kg M**-2 Sec**-1","Kg M**-2 Day**-1",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Pa              ",&
         "Mb              ","Cb              ","Dynes Cm**-2    ",&
         "Mb-1000         ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Ln(Pa)          ","Ln(Mb)          ","Ln(Cb)          ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Pa/Sec          ","Mb/Sec          ",&
         "Mb/Day          ","Cb/Sec          ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Log(Pa)/Sec     ",&
         "Log(Mb)/Sec     ","Log(Cb)/Sec     ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "W M**-2         ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","M**2 Sec**-2    ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Sec/M           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Kg M**-3        ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Kg M**-1 Sec**-1","Kg M**-1 Day**-1",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Kg Sec**-1      ",&
         "10**9 Kg Sec**-1","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           ",&
         "Unset           ","Unset           ","Unset           " /)
    
    !     write out tables
    IF(ifprt(51).GE.2)WRITE(nfprt,1410)
    mm=0
    DO m=1,mxavl
       IF(.NOT.dodia(m))CYCLE
       mm=mm+1
       IF(ifprt(51).GE.2)WRITE(nfprt,1420)mm,m,avail(m)
       IF(ifprt(51).GE.2)WRITE(nfprt,1430)lvavl(m),aunits(nuavl(m)), &
            typcd(itavl(m)),poscd(jpavl(m))
       IF(iavrq(m).NE.0.AND.ifprt(51).GE.2)WRITE(nfprt,1440)iavrq(m)
       IF(kravl(ixavl(m)).LT.0.OR.inavl(m).GT.1)THEN
          ka=ixavl(m)
          kka=ka+inavl(m)-1
          IF(kravl(ka).GT.0)ka=ka+1
          IF(ifprt(51).GE.2)WRITE(nfprt,1450)
          DO kk=ka,kka
             IF(ifprt(51).GE.2)WRITE(nfprt,1460)combf(abs(kravl(kk)))
          END DO
       END IF
    END DO
    IF(ifprt(51).GE.2)WRITE(nfprt,1510)
    DO n=1,mxrq
       IF(ifprt(51).GE.2)WRITE(nfprt,1520)n,reqdg(n)
       IF(ifprt(51).GE.2)WRITE(nfprt,1530)lvrq(n),aunits(nurq(n)), &
            irqcf(n)
       IF(irqav(n).GT.0)THEN
          IF(iclcd(n) < 0 ) THEN
             IF(ifprt(51).GE.2)WRITE(nfprt,1545)irqav(n),reqdg(abs(iclcd(n)))
          ELSE IF (iclcd(n) .EQ. 0 ) THEN
             IF(ifprt(51).GE.2)WRITE(nfprt,1555)irqav(n)
          ELSE
             IF(ifprt(51).GE.2)WRITE(nfprt,1565)irqav(n),reqdg(iclcd(n))
          END IF
       ELSE IF(irqav(n).LT.0)THEN
          irix=abs(irqav(n))
          IF(iclcd(n).LT.0)THEN
             IF(ifprt(51).GE.2)WRITE(nfprt,1575)irix,reqdg(abs(iclcd(n)))
          ELSE
             IF(ifprt(51).GE.2)WRITE(nfprt,1595)irix,reqdg(iclcd(n))
          END IF
       END IF
    END DO
    IF(icf.NE.0)THEN
       IF(ifprt(51).GE.3)WRITE(nfprt,1610)
       DO ix=1,icf
          IF(ifprt(51).GE.3)WRITE(nfprt,1620)ix,combf(ix)
          IF(ifprt(51).GE.3)WRITE(nfprt,1630)lvcf(ix),aunits(nucf(ix)), &
               kfrq(ix),poscd(jpcf(ix))
          IF(ifprt(51).GE.3)WRITE(nfprt,1640)combf(ix)
          ka=ixcf(ix)
          DO kk=ka,incf(ix)+ka-1
             kx=krcf(kk)
             IF(kx.GT.0)ocf=avail(kx)
             IF(kx.LT.0)ocf=combf(-kx)
             jx=jrcf(kk)
             IF(iclcd(jx).LT.0)THEN
                IF(ifprt(51).GE.3)WRITE(nfprt,1650)ocf,kx,jx
             ELSE IF(iclcd(jx).GT.0)THEN
                IF(ifprt(51).GE.3)WRITE(nfprt,1660)ocf,kx,jx
             END IF
          END DO
          IF(icfu(ix))THEN
             IF(ifprt(51).GE.3)WRITE(nfprt,1670)
             ka=ixucf(ix)
             DO kk=ka,inucf(ix)+ka-1
                kx=abs(jrucf(kk))
                IF(ifprt(51).GE.3)WRITE(nfprt,1680)combf(kx)
             END DO
          END IF
       END DO
    END IF

    !     set diagnostic accumulators

    CALL setdia()

    pfbar=iavrq(1).NE.0.OR.iavrq(2).NE.0.OR.iavrq(3).NE.0.OR. &
         iavrq(4).NE.0.OR.iavrq(5).NE.0

30  FORMAT(//'0AVAILABLE DIAGNOSTIC INPUT DECK'/)
50  FORMAT(A40,4I5)
60  FORMAT(' ',A40,4I5)
130 FORMAT(//'0DESIRED DIAGNOSTIC INPUT DECK'/)
150 FORMAT(A40,3I5)
160 FORMAT(' ',A40,3I5)
920 FORMAT(A32)
925 FORMAT(A16)
930 FORMAT(5E16.8)
940 FORMAT(20I4)
1410 FORMAT(' A V A I L A B L E  D I A G N O S T I C S  U S E D  I N  T H I S  R U N')
1420 FORMAT(' DIAG. NO.',I3,' AVAILABLE DIAG. NO.',I3,' NAME = ',A40)
1430 FORMAT(' NUMBER OF LEVELS=',I2,' UNITS: ',A16/' TYPE=',A8,1X,A38)
1440 FORMAT(' REQUESTED DIRECTLY BY DESIRED DIAGNOSTIC NUMBER',I3)
1450 FORMAT(' USED IN THE FOLLOWING COMBINED FIELDS:')
1460 FORMAT(' ',A40)
1510 FORMAT(' D E S I R E D  D I A G N O S T I C S  T A B L E')
1520 FORMAT(' DESIRED DIAGNOSTIC NUMBER',I3,' NAME = ',A40)
1530 FORMAT(' NO. OF LVLS=',I2,' UNITS: ',A16, &
         ' COMBINED FIELD CODE=',I3)
1545 FORMAT(' IS AVAILABLE DIAGNOSTIC NUMBER',I3, &
         ' AND IS SUBTRACTED TO FORM COMBINED FIELD: ',A40)
1555 FORMAT(' IS AVAILABLE DIAGNOSTIC NUMBER',I3,' SAVED DIRECTLY')
1565 FORMAT(' IS AVAILABLE DIAGNOSTIC NUMBER',I3, &
         ' AND IS ADDED TO FORM COMBINED FIELD: ',A40)
1575 FORMAT(' IS COMBINED FIELD NUMBER',I3, &
         ' AND IS SUBTRACTED TO FORM COMBINED FIELD: ',A40)
1595 FORMAT(' IS COMBINED FIELD NUMBER',I3, &
         ' AND IS ADDED TO FORM COMBINED FIELD: ',A40)
1610 FORMAT('1C O M B I N E D  F I E L D  T A B L E')
1620 FORMAT('0COMBINED FIELD NUMBER',I3,' NAME = ',A40)
1630 FORMAT(' NUMBER OF LEVELS=',I2,' UNITS: ',A16/ &
         ' CORRESP. DESIRED DIAG. NO.',I3,1X,A38)
1640 FORMAT(' IS CONSTRUCTED AS FOLLOWS:'/'   ',A40,'=')
1650 FORMAT(' - ',A40,' ( A.D. OR C.F. NO. =',I3,' D.D. NO. =',I3,' )')
1660 FORMAT(' + ',A40,' ( A.D. OR C.F. NO. =',I3,' D.D. NO. =',I3,' )')
1670 FORMAT(' IS USED AS A COMPONENT IN THE FOLLOWING COMBINED FIELDS')
1680 FORMAT(' ',A40)
2100 FORMAT(' DESIRED DIAGNOSTIC DECK NOT WELL ORDERED.'/ &
         ' THE CALCULATION CODE IS ZERO FOR N=',I3, &
         ' WHICH EXCEEDS THE EXPECTED OUTPUT FIELD COUNT =',I3)
2600 FORMAT(' DESIRED DIAGNOSTIC DECK NOT WELL ORDERED.'/ &
         ' THE CALCULATION CODE IS ',I3, &
         ' WHICH EXCEEDS THE CURRENT FIELD COUNT =',I3)
3100 FORMAT(' DESIRED DIAGNOSTIC DECK NOT WELL ORDERED.'/ &
         ' THE CALCULATION CODE IS ',I3,' FOR N=',I3, &
         ' POINTS TO A FIELD WITH NONZERO CALCULATION CODE =',I3)
3600 FORMAT(' DESIRED DIAGNOSTIC DECK NOT WELL ORDERED.',/, &
         ' A FIELD CANNOT BE COMBINED WITH ITSELF.',/, &
         ' THE CALCULATION CODE ',I3,' AND N=',I3, ' POINT TO THE SAME FIELD =',/,'  ',A40)
4100 FORMAT(' DESIRED DIAGNOSTIC FIELD NUMBER',I3,' NAMED ',A40/ &
         ' CANNOT BE FOUND IN AVAILABLE DIAGNOSTICS AND IS NOT REFERENCED AS A COMBINED FIELD')
4600 FORMAT(' DESIRED DIAGNOSTIC FIELD NUMBER',I3,' NAMED ',A40/ &
         ' CANNOT BE FOUND IN AVAILABLE DIAGNOSTICS AND DOES NOT REFERENCE A VALID COMBINED FIELD')
5100 FORMAT(' DESIRED DIAGNOSTIC FIELD NUMBER',I3,' NAMED ',A40/ &
         ' IS NOT USED ANYWHERE')
5600 FORMAT(' COMBINED FIELD COMPONENT TABLE EXCEEDED FOR FIELD NO.', &
         I3,' AND DESIRED DIAGNOSTIC NUMBER',I3)
6100 FORMAT(' AVAILABLE DIAGNOSTIC USE TABLE EXCEEDED FOR NO.', &
         I3,' AND COMBINED FIELD NUMBER',I3)
6600 FORMAT(' AVAILABLE DIAGNOSTIC TABLE EMPTY OR NOT FOUND')
7100 FORMAT(' DESIRED DIAGNOSTIC TABLE EMPTY OR NOT FOUND') 
7600 FORMAT(' TYPE CODES FOR COMBINED FIELD COMPONENTS ARE INCONSISTENT.',/,&
         ' COMBINED FIELD NO.=',I3,' NAME=',A40 / &
         ' AVAILABLE DIAGNOSTIC NO.=',(' ',10I5/))
8100 FORMAT(' COMBINED FIELD USE TABLE EXCEEDED FOR FIELD NO.', &
         I3,' AND COMBINED FIELD NUMBER',I3)
8600 FORMAT(' DESIRED DIAGNOSTIC DECK NOT WELL ORDERED.'/ &
         ' THE COMBINED FIELD',I3,' < COMBINED FIELD COMPONENT',I3/ &
         ' WHICH ARE DESIRED DIAGNOSTICS',I3,' AND',I3)
9100 FORMAT(' UNIT CODE GROUP FOR UNIT CODE',I4, &
         ' OF AVAILABLE DIAGNOSTIC COMPONENT',I3/ &
         ' IS NOT THE SAME CODE GROUP FOR UNIT CODE',I4, &
         ' OF COMBINED FIELD',I3)
9600 FORMAT(' THE UNIT CODE,',I4, &
         ', FOR COMBINED FIELD COMPONENT',I3/ &
         ' IS NOT THE SAME UNIT CODE,',I4, &
         ', FOR COMBINED FIELD',I3)

  END SUBROUTINE InitDiagnostics


  ! pwater :perfoms vertical integration of water vapour.


  SUBROUTINE pwater(q, pwtr, ps, del, nx,mx, nq)
    !
    !
    !  nx......imx=imax+1 or imax+2   :this dimension instead of imax
    !          is used in order to avoid bank conflict of memory
    !          access in fft computation and make it efficient. the
    !          choice of 1 or 2 depends on the number of banks and
    !          the declared type of grid variable (real*4,real*8)
    !          to be fourier transformed.
    !          cyber machine has the symptom.
    !          cray machine has no bank conflict, but the argument
    !          'imx' in subr. fft991 cannot be replaced by imax 
    !  nq......Number of sigma levels      
    !  q.......gq        specific humidity (fourier).
    !          gqu       u*q
    !          gqv       v*q
    !  pwtr
    !  ps......gpphi(imx)            input : latitudinal derivative of
    !                                        natural ig of surface
    !                                        pressure (fourier) 
    !  del.....sigma spacing for each layer computed in routine "setsig". 
    !  grav....grav   gravity constant        (m/s**2)    
    !
    INTEGER, INTENT(IN   ) :: nx
    INTEGER, INTENT(IN   ) :: mx    
    INTEGER, INTENT(IN   ) :: nq   
    REAL(KIND=r8),    INTENT(IN   ) :: q(nx,nq)
    REAL(KIND=r8),    INTENT(INOUT) :: pwtr(nx)
    REAL(KIND=r8),    INTENT(IN   ) :: ps(nx)
    REAL(KIND=r8),    INTENT(IN   ) :: del(nq)  

    REAL(KIND=r8)    :: fac
    INTEGER :: i
    INTEGER :: k

    fac = 1.0e3_r8/grav

    DO i = 1,mx
       pwtr(i) = 0.0_r8
    END DO

    DO k = 1,nq
       DO i = 1,mx
          pwtr(i) = pwtr(i) + q(i,k)*del(k)
       END DO
    END DO

    DO i=1,mx
       pwtr(i) = pwtr(i) * ps(i) * fac
    END DO
  END SUBROUTINE pwater


  SUBROUTINE globme(a, z, dthet, costhe, cf, title, nufr, nuto)
    !
    ! globme :perfoms zonal and global mean.
    !
    !
    !     find global mean of a
    !
    REAL(KIND=r8)             , INTENT(IN   ) :: a(iMaxNew,jMaxNew)
    REAL(KIND=r8)             , INTENT(OUT  ) :: z(jMaxNew)
    REAL(KIND=r8)             , INTENT(IN   ) :: dthet(jMaxNew)
    REAL(KIND=r8)             , INTENT(IN   ) :: costhe(jMaxNew)
    REAL(KIND=r8)             , INTENT(IN   ) :: cf   
    CHARACTER(LEN=40), INTENT(IN   ) :: title
    INTEGER          , INTENT(IN   ) :: nufr   
    INTEGER          , INTENT(IN   ) :: nuto
    INTEGER                          :: i
    INTEGER                          :: j
    REAL(KIND=r8)                             :: gm
    REAL(KIND=r8)  :: work(iMaxNew,jMaxNew)


    DO j=1,jMaxNew
       DO i=1,iMaxNew
          work(i,j)=a(i,j)*cf
       END DO
    END DO
    CALL cnvray(work  ,iMaxNew*jMaxNew, nufr, nuto)
    !
    !     zonal mean
    !
    DO j=1,jMaxNew
       z(j) = 0.0_r8
       DO i=1,iMaxNew
          z(j) = z(j) + work(i,j)
       END DO
       z(j) = z(j)/float(iMaxNew)
    END DO
    !
    !     integral with latitude
    !
    gm = 0.0_r8
    DO j=1,jMaxNew
       gm = gm + z(j)*costhe(j)*dthet(j)
    END DO
    gm = 0.5_r8 * gm
    IF(ifprt(31).GE.1)WRITE(nfprt,70)title,aunits(nuto),gm
    !
    !     global mean standard deviation
    !
    DO j=1,jMaxNew
       DO i=1,iMaxNew
          work(i,j)=(work(i,j)-gm)*(work(i,j)-gm)
       END DO
    END DO

    DO j=1,jMaxNew
       z(j) = 0.0_r8
       DO i=1,iMaxNew
          z(j) = z(j) + work(i,j)
       END DO
       z(j) = z(j)/float(iMaxNew)
    END DO
    !
    !     integral with latitude
    !
    gm = 0.0_r8
    DO j=1,jMaxNew
       gm = gm + z(j)*costhe(j)*dthet(j)
    END DO
    gm = SQRT(0.5_r8 * gm)
    WRITE(6,71)title,aunits(nuto),gm
70  FORMAT(' ',A40,' IN UNITS OF ',A16,1X,G16.8,' G.M.')
71  FORMAT(' ',A40,' IN UNITS OF ',A16,1X,G16.8,' S.D.')
  END SUBROUTINE globme












  SUBROUTINE setdia()
    !
    ! setdia :extended diagnostics version 1 used for
    !         initializing and partitioning the diagnostic accumulators;
    !         this version is the memory resident version;
    !         see subroutine indiag for further discussion.
    !
    INTEGER :: nf
    INTEGER :: ix
    INTEGER :: m 
    INTEGER :: nw 
    INTEGER :: ish
    INTEGER :: igh

    mxgaus=mgaus*kMaxNew+ngaus
    mxspec=mspec*kMaxNew+nspec

    ALLOCATE (spec1d(ibMax, mxgaus, jbMax))
    spec1d=0.0_r8
    
    ALLOCATE (spec(2*mnMax, mxspec))
    spec = 0.0_r8

    ALLOCATE (gaus(ibMax, mxgaus, jbMax))
    gaus = 0.0_r8

    ALLOCATE (gaus_in(iMaxNew, mxgaus, jMaxNew))
    gaus_in = 0.0_r8

    ALLOCATE (CountGaus(ibMax, mxgaus, jbMax))
    CountGaus = 0.0_r8

    ALLOCATE (CountTOTAL(mxgaus, jbMax))
    CountTOTAL = 0.0_r8

    ispec=0
    igaus=0
    nf=0
    ix=0
    DO m=1,mxavl
       lspec(m)=0
       lgaus(m)=0
       IF(.NOT.dodia(m))THEN
          CONTINUE
       ELSE IF(iavrq(m).LE.0)THEN
          CONTINUE
       ELSE IF(itavl(m).EQ.1)THEN
          nf=nf+1
          lgaus(m)=igaus+1
          igaus=igaus+lvavl(m)
          IF(igaus.GT.mxgaus)THEN
             WRITE(nfprt,20100)igaus,m,ix
             WRITE(nferr,20100)igaus,m,ix
             STOP 20100
          END IF
       ELSE IF(itavl(m).EQ.2)THEN
          nf=nf+1
          lspec(m)=ispec+1
          ispec=ispec+lvavl(m)
          IF(ispec.GT.mxspec)THEN
             WRITE(nfprt,20600)ispec,m,ix
             WRITE(nferr,20600)ispec,m,ix
             STOP 20600
          END IF
       END IF
    END DO

    IF(nf+icf.NE.nof)THEN
       WRITE(nfprt,21100)nf,icf,nof
       WRITE(nferr,21100)nf,icf,nof
       STOP 21100
    END IF

    IF(ifprt(71).GE.1)WRITE(nfprt,200)nf,ispec,igaus
    IF(icf.NE.0)THEN
       ish=ispec
       igh=igaus
       DO ix=-1,-icf,-1
          IF(itcf(-ix).EQ.1)THEN
             nf=nf+1
             lgaus(ix)=igaus+1
             igaus=igaus+lvcf(-ix)
             IF(igaus.GT.mxgaus)THEN
                WRITE(nfprt,20100)igaus,m,ix
                WRITE(nferr,20100)igaus,m,ix
                STOP 20100
             END IF
          ELSE IF(itcf(-ix).EQ.2)THEN
             nf=nf+1
             lspec(ix)=ispec+1
             ispec=ispec+lvcf(-ix)
             IF(ispec.GT.mxspec)THEN
                WRITE(nfprt,20600)ispec,m,ix
                WRITE(nferr,20600)ispec,m,ix
                STOP 20600
             END IF
          END IF
       END DO
       IF(ifprt(71).GE.1)WRITE(nfprt,400)icf,ispec-ish,igaus-igh
    END IF
    nw=ispec*2*mnMax+igaus*ibMax*jbMax
    IF(ifprt(71).GE.1)WRITE(nfprt,500)nf,ispec,igaus,nw

200 FORMAT(' ACCUMULATORS FOR DIRECT DIAGNOSTICS ARE:'/' ',I3, &
         ' FIELDS USING',I4,' SPECTRAL SLOTS AND',I4, &
         ' GAUSSIAN SLOTS')
400 FORMAT('0ACCUMULATORS FOR COMBINED FIELDS ARE:'/' ',I3, &
         ' FIELDS USING',I4,' SPECTRAL SLOTS AND',I4, &
         ' GAUSSIAN SLOTS')
500 FORMAT(' TOTAL ACCUMULATORS FOR ALL DIAGNOSTICS ARE:'/' ',I3, &
         ' FIELDS USING',I4,' SPECTRAL SLOTS AND',I4,' GAUSSIAN SLOTS' &
         ,/,' TOTAL WORDS IN DIAGNOSTIC ACCUMULATORS=',I16//)
20600 FORMAT(' TOTAL NUMBER OF AVAILABLE SPECTRAL SLOTS EXCEEDED.'/ &
         ' NO. OF SLOTS AT LIMIT POINT:',I4,' M=',I4,' IX=',I4)
20100 FORMAT(' TOTAL NUMBER OF AVAILABLE GAUSSIAN SLOTS EXCEEDED.'/ &
         ' NO. OF SLOTS AT LIMIT POINT:',I4,' M=',I4,' IX=',I4)
21100 FORMAT(' NUMBER OF COMPUTED FIELDS =',I3,' +',I3, &
         ' IS NOT THE SAME AS THE NUMBER OF DESIRED FIELDS=',I4)
  END SUBROUTINE setdia





  SUBROUTINE rsdiag()
    !
    ! rsdiag :extended diagnostics version 1;
    !         reset all the diagnostic accumulators;
    !         this version is the memory resident version;
    !         see subroutine indiag for further discussion.
    !
    IF (ispec > 0) THEN
       spec = 0.0_r8
    END IF
    IF (igaus > 0) THEN
       gaus = 0.0_r8
    END IF
  END SUBROUTINE rsdiag

  SUBROUTINE getgau(field, jj)
    !
    !
    ! getgau :extended diagnostics version 1 diagnostic field fetching routine
    !         memory resident version;  see subroutine indiag for further
    !         discussion; used for gaussian fields only
    !
    !
    !     version 1 diagnostic accumulators
    !
    REAL(KIND=r8)   , INTENT(OUT) :: field(iMaxNew,jMaxNew,kMaxNew)
    INTEGER, INTENT(IN ) :: jj 

    INTEGER :: kg  
    INTEGER :: lvl
    INTEGER :: j
    INTEGER :: i
    INTEGER :: l 
    INTEGER :: ll

    IF (jj == 0) THEN
       WRITE(nfprt,3270)
       WRITE(nferr,3270)
       STOP 3270
    ELSE IF (jj >  0) THEN
       lvl=lvavl(jj)
    ELSE IF (jj <  0) THEN
       lvl=lvcf(-jj)
    END IF
    kg=lgaus(jj)
    DO j = 1, jMaxNew
       DO l = 1, lvl
          ll = kg+l-1
          DO i = 1, iMaxNew
             field(i,j,l)=gaus_in(i,ll,j)
          END DO
       END DO
    END DO
3270 FORMAT(' ERROR IN CALLING GETGAU WITH NO SET LOCATION')
  END SUBROUTINE getgau






  SUBROUTINE wridia (iudiag, maxstp)
    !
    ! wridia :writes kistler/katz/schneider diagnostic fields on disk.
    !
    INTEGER, INTENT(IN) :: iudiag
    INTEGER, INTENT(IN) :: maxstp

    REAL(KIND=r8)    :: gwork(iMaxNew,jMaxNew,kMaxNew)
    REAL(KIND=r8)    :: zonmn(jMaxNew)
    LOGICAL(KIND=i8) :: fgm
    REAL(KIND=r8)    :: f1
    INTEGER :: j 
    INTEGER :: m 
    INTEGER :: nn
    INTEGER :: mm
    INTEGER :: ii
    INTEGER :: ix
    INTEGER :: ka
    INTEGER :: ki
    INTEGER :: kk
    INTEGER :: kx
    INTEGER :: jx
    INTEGER :: ji
    INTEGER :: l 
    INTEGER :: i
    INTEGER :: k     
    !     
    !     global mean diagnostics printed when available
    !     
    !     name:                                       a.d. no.
    !     total precipitation                             9
    !     convective precipitation                       10
    !     surface sensible heat flux                     16
    !     surface zonal wind stress                      18
    !     surface meridional wind stress                 19
    !     surface latent heating                         17
    !     downward longwave flux at the bottom           21
    !     outgoing longwave radiation at the top         23
    !     incident shortwave flux at the top             24
    !     downward shorwave at the ground                25
    !     shortwave absorbed by the earth/atmosphere     28
    !     shortwave absorbed at the ground               29
    !     net longwave at the ground                     30
    !     
    !
    !     f1 for time intensive (1/maxstp)
    !     synoptic interval is 24 hours presumeably but can be any integral
    !     number of time steps
    !
    f1=1.0_r8/float(maxstp)
    fgm=.TRUE.
    gaus_in=gaus
    !     
    !     directly available fields
    !     
    DO m=1,mxavl
       IF ((dodia(m)) .AND. (iavrq(m) > 0)) THEN
          nn=iavrq(m)
          IF (itavl(m) == 2) THEN
             mm=lspec(m)
             !CALL sclout(iudiag, spec(1,mm), ihdim(itavl(m)), lvavl(m), f1, &
             !    ivar(itavl(m)), nuavl(m), nurq(nn))
          ELSE
             CALL getgau(gwork, m)
             CALL sclout(iudiag, gwork, ihdim(itavl(m)), lvavl(m), &
                  f1, ivar(itavl(m)), nuavl(m), nurq(nn))
             DO ii=1,ngbme
                IF (igbme(ii) == m) GO TO 250
             END DO
             CYCLE
250          IF (fgm .AND. ifprt(31) >= 1)WRITE(nfprt,300)
             IF (fgm .AND. dodia(17) .AND. ifprt(31) >= 2)WRITE(nfprt,350)
             CALL globme(gwork, zonmn, dcol, scol, f1, avail(m), &
                  nuavl(m), nurq(nn))
             fgm=.FALSE.
          END IF
       END IF
    END DO
    IF (.NOT. fgm .AND. ifprt(31) >= 1) WRITE(nfprt,1010) 
    !
    !     combined fields
    !
    IF (icf /= 0) THEN
       !
       !     first combine other combined field components
       !
       DO ix=1,icf
          IF (icfu(ix)) THEN
             ka=ixucf(ix)
             IF (itcf(ix) == 1) THEN
                ki=lgaus(-ix)
                DO j=1,jMaxNew
                   DO kk=ka,inucf(ix)+ka-1
                      kx=jrucf(kk)
                      jx=abs(kx)
                      ji=lgaus(-jx)
                      DO l=1,lvcf(jx)
                         IF (kx > 0) THEN
                            DO i=1,iMaxNew
                               gaus_in(i,l+ji-1,j)=gaus_in(i,l+ji-1,j)+ &
                                    gaus_in(i,l+ki-1,j)
                            END DO
                         ELSE
                            DO i=1,iMaxNew
                               gaus_in(i,l+ji-1,j)=gaus_in(i,l+ji-1,j)- &
                                    gaus_in(i,l+ki-1,j)
                            END DO
                         END IF
                      END DO
                   END DO
                END DO
             ELSE IF (itcf(ix) == 2) THEN
                ki=lspec(-ix)
                DO kk=ka,inucf(ix)+ka-1
                   kx=jrucf(kk)
                   jx=abs(kx)
                   ji=lspec(-jx)
                   DO l=1,lvcf(jx)
                      IF (kx > 0) THEN
                         DO i=1,2*mnMax
                            spec(i,l+ji-1)=spec(i,l+ji-1)+ &
                                 spec(i,l+ki-1)
                         END DO
                      ELSE
                         DO i=1,2*mnMax
                            spec(i,l+ji-1)=spec(i,l+ji-1)- &
                                 spec(i,l+ki-1)
                         END DO
                      END IF
                   END DO
                END DO
             END IF
          END IF
       END DO
       !
       !     obtain combined fields and write out
       !
       DO ix=1,icf
          IF (itcf(ix) == 2) THEN
             mm=lspec(-ix)
!             CALL sclout(iudiag,spec(1,mm),ihdim(itcf(ix)),lvcf(ix), &
!                  f1,ivar(itcf(ix)),nucf(ix),nucf(ix))
          ELSE IF (itcf(ix) == 1)THEN
             CALL getgau(gwork, -ix)
              WRITE(*,*)avail(-ix)
             CALL sclout(iudiag,gwork,ihdim(itcf(ix)),lvcf(ix), &
                  f1,ivar(itcf(ix)),nucf(ix),nucf(ix))
          END IF
       END DO
    END IF
    gaus_in=0.0
300 FORMAT(//'0GLOBAL MEAN DIAGNOSTICS'//)
350 FORMAT(' NOTE: TO COMPUTE EVAPORATION IN MM/DAY FROM LATENT '  &
         ,'HEATING DIVIDE BY 28.9'//)
1010 FORMAT(//)
  END SUBROUTINE wridia






  SUBROUTINE updiaA(field, loca, lat)
    !
    ! updia  :extended diagnostics version 1 diagnostic field accumulator
    !         subroutine; memory resident version;
    !         see subroutine indiag for further discussion;
    !         for gaussian fields only called one gaussian latitude at a time.
    !
    IMPLICIT NONE
    REAL(KIND=r8),    INTENT(in   ) :: field(:)
    INTEGER, INTENT(in   ) :: loca
    INTEGER, INTENT(in   ) :: lat
    REAL(KIND=r8)     :: hold(ibMax)
    INTEGER  :: imkm 
    INTEGER  :: i
    INTEGER  :: l 
    INTEGER  :: ll
    INTEGER  :: ka 
    INTEGER  :: kg 
    INTEGER  :: kka 
    INTEGER  :: kk 
    INTEGER  :: kcf 
    INTEGER  :: kx
    INTEGER  :: ja 
    INTEGER  :: jj 
    INTEGER  :: jja 
    INTEGER  :: jx 
    INTEGER  :: jcf

    IF (.NOT. dodia(loca)) THEN
       WRITE(nfprt,3180)loca
       WRITE(nferr,3180)loca
       STOP 3180
    END IF

    IF (itavl(loca) /= 1) THEN
       WRITE(nfprt,4180)itavl(loca)
       WRITE(nferr,4180)itavl(loca)
       STOP 4180
    END IF
    ka=ixavl(loca)
    imkm = ibMax
    !
    !    case for directly saved fields
    !
    IF (iavrq(loca) > 0) THEN
       kg=lgaus(loca)
       DO i = 1, ibMaxPerJB(lat)
          gaus(i,kg,lat)=gaus(i,kg,lat)+field(i)
       END DO
    END IF
    !
    !    case for combined fields
    !
    IF (kravl(ka) < 0 .OR. inavl(loca) > 1) THEN
       kka=ka+inavl(loca)-1
       IF (kravl(ka) > 1) ka=ka+1
       !
       !    for each combined field using the supplied available diagnostic
       !    
       DO kk = ka, kka
          kcf=kravl(kk)
          jcf=-kcf
          ja =ixcf(jcf)
          jja=ja+incf(jcf)-1
          !    
          !    search for corresponding desired field
          !    
          DO jj = ja, jja
             jx=jrcf(jj)
             kx=krcf(jj)
             IF (kx == loca) go to 200
          END DO

          WRITE(nfprt,3680)loca,jcf,kk,ja,jja
          WRITE(nferr,3680)loca,jcf,kk,ja,jja

          STOP 3680
          !    
          !    treat each accumulation according the the sign of the desired
          !    calculation code (iclcd)
          !    
200       CONTINUE
          DO i = 1, ibMaxPerJB(lat)
             hold(i)=field(i)
          END DO

          CALL cnvray(hold,imkm,nuavl(loca),nucf(jcf))

          IF (iclcd(jx) < 0) THEN
             kg=lgaus(kcf)
             DO i = 1, 1
                gaus(i,kg,lat)=gaus(i,kg,lat)-hold(i)
             END DO
          ELSE
             kg=lgaus(kcf)
             DO i = 1, 1
                gaus(i,kg,lat)=gaus(i,kg,lat)+hold(i)
             END DO
          END IF
       END DO
    END IF
3180 FORMAT(' ERROR IN CALLING UPDIA WITH UNSET AVAILABLE DIAGNOSTIC', I3)
3680 FORMAT(' UNABLE TO FIND MATCHING AVAILABLE DIAG. NO.',I3/ &
         ' FOR COMBINED FIELD',I3,' A.D. INDEX',I3,' C.F. RANGE',I3,'-',I3)
4555 FORMAT(' CONVERSION ERROR IN UPDIA.  ERROR=',I3,' NUAVL=',I5, &
         ' NUCF=',I5/' A.D. NO.=',I3,' C.F. NO.=',I3,' A.D. INDEX=',I3)
4180 FORMAT(' ERROR IN CALLING UPDIA WITH WRONG TYPE CODE',I2)
  END SUBROUTINE updiaA



  SUBROUTINE StoreMaskedDiag1D(field, loca, lat,jdt)
    !
    ! updia  :extended diagnostics version 1 diagnostic field accumulator
    !         subroutine; memory resident version;
    !         see subroutine indiag for further discussion;
    !         for gaussian fields only called one gaussian latitude at a time.
    !
    IMPLICIT NONE
    REAL(KIND=r8),    INTENT(in   ) :: field(:)
    INTEGER, INTENT(in   ) :: loca
    INTEGER, INTENT(in   ) :: lat
    INTEGER, INTENT(in   ) :: jdt
    REAL(KIND=r8)     :: hold(ibMax)
    INTEGER  :: imkm 
    INTEGER  :: i
    INTEGER  :: l 
    INTEGER  :: ll
    INTEGER  :: ka 
    INTEGER  :: kg 
    INTEGER  :: kka 
    INTEGER  :: kk 
    INTEGER  :: kcf 
    INTEGER  :: kx
    INTEGER  :: ja 
    INTEGER  :: jj 
    INTEGER  :: jja 
    INTEGER  :: jx 
    INTEGER  :: jcf

    IF (.NOT. dodia(loca)) THEN
       WRITE(nfprt,3180)loca
       WRITE(nferr,3180)loca
       STOP 3180
    END IF

    IF (itavl(loca) /= 1) THEN
       WRITE(nfprt,4180)itavl(loca)
       WRITE(nferr,4180)itavl(loca)
       STOP 4180
    END IF
    ka=ixavl(loca)
    imkm = ibMax
    !
    !    case for directly saved fields
    !
    IF (iavrq(loca) > 0) THEN
       kg=lgaus(loca)
       DO i = 1, ibMaxPerJB(lat)
          CountTOTAL(kg,lat)=CountTOTAL(kg,lat)+1.0_r8   !nilo
          IF (field(i) /= undef) THEN
             gaus     (i,kg,lat)=gaus     (i,kg,lat)+field(i)
             CountGaus(i,kg,lat)=CountGaus(i,kg,lat)+1.0_r8
         END IF
         IF(cthl(jdt))THEN
            IF (CountGaus(i,kg,lat) /= 0.0_r8 ) THEN
                gaus     (i,kg,lat) = CountTOTAL(kg,lat)*(gaus(i,kg,lat)/CountGaus(i,kg,lat))
            ELSE
                gaus     (i,kg,lat) = 0.0_r8
            END IF
         END IF
       END DO
    END IF
    !
    !    case for combined fields
    !
    IF (kravl(ka) < 0 .OR. inavl(loca) > 1) THEN
       kka=ka+inavl(loca)-1
       IF (kravl(ka) > 1) ka=ka+1
       !
       !    for each combined field using the supplied available diagnostic
       !    
       DO kk = ka, kka
          kcf=kravl(kk)
          jcf=-kcf
          ja =ixcf(jcf)
          jja=ja+incf(jcf)-1
          !    
          !    search for corresponding desired field
          !    
          DO jj = ja, jja
             jx=jrcf(jj)
             kx=krcf(jj)
             IF (kx == loca) go to 200
          END DO

          WRITE(nfprt,3680)loca,jcf,kk,ja,jja
          WRITE(nferr,3680)loca,jcf,kk,ja,jja

          STOP 3680
          !    
          !    treat each accumulation according the the sign of the desired
          !    calculation code (iclcd)
          !    
200       CONTINUE
          DO i = 1, ibMaxPerJB(lat)
             hold(i)=field(i)
          END DO

          CALL cnvray(hold,imkm,nuavl(loca),nucf(jcf))

          IF (iclcd(jx) < 0) THEN
             kg=lgaus(kcf)
             DO i = 1, ibMaxPerJB(lat)
                IF (field(i) /= undef) THEN
                   gaus     (i,kg,lat)=gaus     (i,kg,lat)-hold(i)
                   CountGaus(i,kg,lat)=CountGaus(i,kg,lat)+1.0_r8
                ELSE
                   gaus     (i,kg,lat)=gaus     (i,kg,lat)
                   CountGaus(i,kg,lat)=CountGaus(i,kg,lat)
                END IF
                IF(cthl(jdt))THEN
                   IF (CountGaus(i,kg,lat) /= 0.0_r8 ) THEN
                       gaus     (i,kg,lat) = CountTOTAL(kg,lat)*(gaus(i,kg,lat)/CountGaus(i,kg,lat))
                   ELSE
                       gaus     (i,kg,lat) = CountTOTAL(kg,lat)*0.0_r8
                   END IF
                END IF
             END DO
          ELSE
             kg=lgaus(kcf)
             DO i = 1, ibMaxPerJB(lat)
                IF (field(i) /= undef) THEN
                   gaus     (i,kg,lat)=gaus     (i,kg,lat)+hold(i)
                   CountGaus(i,kg,lat)=CountGaus(i,kg,lat)+1.0_r8
                ELSE
                   gaus     (i,kg,lat)=gaus     (i,kg,lat)
                   CountGaus(i,kg,lat)=CountGaus(i,kg,lat)
                END IF
                IF(cthl(jdt))THEN
                   IF (CountGaus(i,kg,lat) /= 0.0_r8 ) THEN
                       gaus     (i,kg,lat) = CountTOTAL(kg,lat)*(gaus(i,kg,lat)/CountGaus(i,kg,lat))
                   ELSE
                       gaus     (i,kg,lat) = CountTOTAL(kg,lat)*0.0_r8
                   END IF
                END IF
             END DO
          END IF
       END DO
    END IF
3180 FORMAT(' ERROR IN CALLING UPDIA WITH UNSET AVAILABLE DIAGNOSTIC', I3)
3680 FORMAT(' UNABLE TO FIND MATCHING AVAILABLE DIAG. NO.',I3/ &
         ' FOR COMBINED FIELD',I3,' A.D. INDEX',I3,' C.F. RANGE',I3,'-',I3)
4555 FORMAT(' CONVERSION ERROR IN UPDIA.  ERROR=',I3,' NUAVL=',I5, &
         ' NUCF=',I5/' A.D. NO.=',I3,' C.F. NO.=',I3,' A.D. INDEX=',I3)
4180 FORMAT(' ERROR IN CALLING UPDIA WITH WRONG TYPE CODE',I2)
  END SUBROUTINE StoreMaskedDiag1D


  SUBROUTINE updiaB(field, loca, lat)
    !
    ! updia  :extended diagnostics version 1 diagnostic field accumulator
    !         subroutine; memory resident version;
    !         see subroutine indiag for further discussion;
    !         for gaussian fields only called one gaussian latitude at a time.
    !
    IMPLICIT NONE
    REAL(KIND=r8),    INTENT(in   ) :: field(:,:)
    INTEGER, INTENT(in   ) :: loca
    INTEGER, INTENT(in   ) :: lat
    REAL(KIND=r8)     :: hold(ibMax,kMaxNew)
    INTEGER  :: imkm 
    INTEGER  :: i
    INTEGER  :: lvl 
    INTEGER  :: l 
    INTEGER  :: ll
    INTEGER  :: ka 
    INTEGER  :: kg 
    INTEGER  :: kka 
    INTEGER  :: kk 
    INTEGER  :: kcf 
    INTEGER  :: kx
    INTEGER  :: ja 
    INTEGER  :: jj 
    INTEGER  :: jja 
    INTEGER  :: jx 
    INTEGER  :: jcf

    IF (.NOT. dodia(loca)) THEN
       WRITE(nfprt,3180)loca
       WRITE(nferr,3180)loca
       STOP 3180
    END IF

    IF (itavl(loca) /= 1) THEN
       WRITE(nfprt,4180)itavl(loca)
       WRITE(nferr,4180)itavl(loca)
       STOP 4180
    END IF
    lvl=lvavl(loca)
    ka=ixavl(loca)
    imkm = ibMax*lvl
    !
    !    case for directly saved fields
    !
    IF (iavrq(loca) > 0) THEN
       kg=lgaus(loca)
       DO l = 1, lvl
          ll=kg+l-1
          DO i = 1, ibMaxPerJB(lat)
             gaus(i,ll,lat)=gaus(i,ll,lat)+field(i,l)
          END DO
       END DO
    END IF
    !
    !    case for combined fields
    !
    IF (kravl(ka) < 0 .OR. inavl(loca) > 1) THEN
       kka=ka+inavl(loca)-1
       IF (kravl(ka) > 1) ka=ka+1
       !
       !    for each combined field using the supplied available diagnostic
       !    
       DO kk = ka, kka
          kcf=kravl(kk)
          jcf=-kcf
          ja =ixcf(jcf)
          jja=ja+incf(jcf)-1
          !    
          !    search for corresponding desired field
          !    
          DO jj = ja, jja
             jx=jrcf(jj)
             kx=krcf(jj)
             IF (kx == loca) go to 200
          END DO

          WRITE(nfprt,3680)loca,jcf,kk,ja,jja
          WRITE(nferr,3680)loca,jcf,kk,ja,jja

          STOP 3680
          !    
          !    treat each accumulation according the the sign of the desired
          !    calculation code (iclcd)
          !    
200       CONTINUE

          DO l = 1, lvl
             DO i = 1, 1
                hold(i,l)=field(i,l)
             END DO
          END DO

          CALL cnvray(hold,imkm,nuavl(loca),nucf(jcf))

          IF (iclcd(jx) < 0) THEN
             kg=lgaus(kcf)
             DO l = 1, lvl
                ll=kg+l-1
                DO i = 1,  ibMaxPerJB(lat)
                   gaus(i,ll,lat)=gaus(i,ll,lat)-hold(i,l)
                END DO
             END DO
          ELSE
             kg=lgaus(kcf)
             DO l = 1, lvl
                ll=kg+l-1
                DO i = 1,  ibMaxPerJB(lat)
                   gaus(i,ll,lat)=gaus(i,ll,lat)+hold(i,l)
                END DO
             END DO
          END IF
       END DO
    END IF
3180 FORMAT(' ERROR IN CALLING UPDIA WITH UNSET AVAILABLE DIAGNOSTIC', I3)
3680 FORMAT(' UNABLE TO FIND MATCHING AVAILABLE DIAG. NO.',I3/ &
         ' FOR COMBINED FIELD',I3,' A.D. INDEX',I3,' C.F. RANGE',I3,'-',I3)
4555 FORMAT(' CONVERSION ERROR IN UPDIA.  ERROR=',I3,' NUAVL=',I5, &
         ' NUCF=',I5/' A.D. NO.=',I3,' C.F. NO.=',I3,' A.D. INDEX=',I3)
4180 FORMAT(' ERROR IN CALLING UPDIA WITH WRONG TYPE CODE',I2)
  END SUBROUTINE updiaB



  SUBROUTINE StoreMaskedDiag2D(field, loca, lat ,jdt)
    !
    ! updia  :extended diagnostics version 1 diagnostic field accumulator
    !         subroutine; memory resident version;
    !         see subroutine indiag for further discussion;
    !         for gaussian fields only called one gaussian latitude at a time.
    !
    IMPLICIT NONE
    REAL(KIND=r8),    INTENT(in   ) :: field(:,:)
    INTEGER, INTENT(in   ) :: loca
    INTEGER, INTENT(in   ) :: lat
    INTEGER, INTENT(in   ) :: jdt
    REAL(KIND=r8)     :: hold(ibMax,kMaxNew)
    INTEGER  :: imkm 
    INTEGER  :: i
    INTEGER  :: lvl 
    INTEGER  :: l 
    INTEGER  :: ll
    INTEGER  :: ka 
    INTEGER  :: kg 
    INTEGER  :: kka 
    INTEGER  :: kk 
    INTEGER  :: kcf 
    INTEGER  :: kx
    INTEGER  :: ja 
    INTEGER  :: jj 
    INTEGER  :: jja 
    INTEGER  :: jx 
    INTEGER  :: jcf
    IF (.NOT. dodia(loca)) THEN
       WRITE(nfprt,3180)loca
       WRITE(nferr,3180)loca
       STOP 3180
    END IF
    IF (itavl(loca) /= 1) THEN
       WRITE(nfprt,4180)itavl(loca)
       WRITE(nferr,4180)itavl(loca)
       STOP 4180
    END IF
    lvl=lvavl(loca)
    ka=ixavl(loca)
    imkm = ibMax*lvl
    !
    !    case for directly saved fields
    !
    IF (iavrq(loca) > 0) THEN
       kg=lgaus(loca)
       DO l = 1, lvl
          ll=kg+l-1
          CountTOTAL(ll,lat)=CountTOTAL(ll,lat)+1.0_r8 !nilo
          DO i = 1, ibMaxPerJB(lat)
             IF (field(i,l) /= undef) THEN
                gaus(i,ll,lat)=gaus(i,ll,lat)+field(i,l)
                CountGaus(i,ll,lat)=CountGaus(i,ll,lat)+1.0_r8
             END IF
             IF (cthl(jdt)) THEN
                IF (CountGaus(i,ll,lat) /= 0.0_r8) THEN
                   gaus(i,ll,lat) = CountTOTAL(ll,lat)*&
                        (gaus(i,ll,lat)/CountGaus(i,ll,lat))
                ELSE
                   gaus(i,ll,lat) = 0.0_r8
                END IF
             END IF
          END DO
       END DO
    END IF
    !
    !    case for combined fields
    !
    IF (kravl(ka) < 0 .OR. inavl(loca) > 1) THEN
       kka=ka+inavl(loca)-1
       IF (kravl(ka) > 1) ka=ka+1
       !
       !    for each combined field using the supplied available diagnostic
       !    
       DO kk = ka, kka
          kcf=kravl(kk)
          jcf=-kcf
          ja =ixcf(jcf)
          jja=ja+incf(jcf)-1
          !    
          !    search for corresponding desired field
          !    
          DO jj = ja, jja
             jx=jrcf(jj)
             kx=krcf(jj)
             IF (kx == loca) go to 200
          END DO

          WRITE(nfprt,3680)loca,jcf,kk,ja,jja
          WRITE(nferr,3680)loca,jcf,kk,ja,jja

          STOP 3680
          !    
          !    treat each accumulation according the the sign of the desired
          !    calculation code (iclcd)
          !    
200       CONTINUE

          DO l = 1, lvl
             DO i = 1, 1
                hold(i,l)=field(i,l)
             END DO
          END DO

          CALL cnvray(hold,imkm,nuavl(loca),nucf(jcf))

          IF (iclcd(jx) < 0) THEN
             kg=lgaus(kcf)
             DO l = 1, lvl
                ll=kg+l-1
                DO i = 1,  ibMaxPerJB(lat)
                   IF (field(i) /= undef) THEN
                      gaus     (i,ll,lat)=gaus     (i,ll,lat)-hold(i,l)
                      CountGaus(i,ll,lat)=CountGaus(i,ll,lat)-1.0_r8
                   ELSE
                      gaus     (i,ll,lat)=gaus     (i,ll,lat)
                      CountGaus(i,ll,lat)=CountGaus(i,ll,lat)
                   END IF
                   IF(cthl(jdt))THEN
                      IF (CountGaus(i,ll,lat) /= 0.0_r8 ) THEN
                          gaus     (i,ll,lat) = CountTOTAL(ll,lat)*(gaus(i,ll,lat)/CountGaus(i,ll,lat))
                      ELSE
                          gaus     (i,ll,lat) = CountTOTAL(ll,lat)*0.0_r8
                     END IF
                   END IF
                END DO
             END DO
          ELSE
             kg=lgaus(kcf)
             DO l = 1, lvl
                ll=kg+l-1
                DO i = 1,  ibMaxPerJB(lat)
                   IF (field(i) /= undef) THEN
                      gaus     (i,ll,lat)=gaus     (i,ll,lat)+hold(i,l)
                      CountGaus(i,ll,lat)=CountGaus(i,ll,lat)+1.0_r8
                   ELSE
                      gaus     (i,ll,lat)=gaus     (i,ll,lat)
                      CountGaus(i,ll,lat)=CountGaus(i,ll,lat)
                   END IF
                   IF(cthl(jdt))THEN
                      IF (CountGaus(i,kg,lat) /= 0.0_r8 ) THEN
                         gaus     (i,ll,lat) = CountTOTAL(ll,lat)*(gaus(i,ll,lat)/CountGaus(i,ll,lat))
                      ELSE
                         gaus     (i,ll,lat) = CountTOTAL(ll,lat)*0.0_r8
                      END IF
                   END IF
                END DO
             END DO
          END IF
       END DO
    END IF
3180 FORMAT(' ERROR IN CALLING UPDIA WITH UNSET AVAILABLE DIAGNOSTIC', I3)
3680 FORMAT(' UNABLE TO FIND MATCHING AVAILABLE DIAG. NO.',I3/ &
         ' FOR COMBINED FIELD',I3,' A.D. INDEX',I3,' C.F. RANGE',I3,'-',I3)
4555 FORMAT(' CONVERSION ERROR IN UPDIA.  ERROR=',I3,' NUAVL=',I5, &
         ' NUCF=',I5/' A.D. NO.=',I3,' C.F. NO.=',I3,' A.D. INDEX=',I3)
4180 FORMAT(' ERROR IN CALLING UPDIA WITH WRONG TYPE CODE',I2)
  END SUBROUTINE StoreMaskedDiag2D




  SUBROUTINE reord (datum, dim2, work, lev, imask, tsea, ittl)
    INTEGER, INTENT(IN ) :: dim2
    REAL(KIND=r8),    INTENT(in ) :: datum(ibMax,dim2,jbMax)
    REAL(KIND=r8),    INTENT(OUT) :: work (ibMax,jbMax)
    INTEGER, INTENT(IN ) :: lev
    INTEGER(KIND=i8), INTENT(IN ) :: imask   (ibMax,jbMax)
    REAL(KIND=r8),    INTENT(IN ) :: tsea    (ibMax,jbMax)
    CHARACTER(LEN=4), INTENT(IN) :: ittl

    INTEGER :: j
    INTEGER :: i
    INTEGER :: ncount
    LOGICAL(KIND=i8) :: case1
    LOGICAL(KIND=i8) :: case2
    case1 = ittl == 'TD  '
    case2 = ittl == 'W1  ' .OR. ittl == 'W2  ' .OR. ittl == 'W3  ' 
    DO j = 1, jbMax
       ncount=0
       DO i = 1, ibMax
          IF (imask(i,j) >= 1) THEN
             ncount = ncount + 1
             work(i,j) = datum(i,lev,j)
          ELSE IF (case1) THEN
             work(i,j)=ABS(tsea(i,j))
          ELSE IF (case2) THEN
             work(i,j)=1.0_r8
          ELSE
             work(i,j)=0.0_r8
          END IF
       END DO
    END DO
  END SUBROUTINE reord






  SUBROUTINE wrfct(nfdrct, nfdiag, nffcst, ifday, tod, idate, idatec,&
                    roperm,namef,labeli,labelf,extw,exdw,trunc,lev,longitude,latitude,&
                    pmand,it)
    INTEGER, INTENT(IN ) :: nfdrct
    INTEGER, INTENT(IN ) :: nfdiag
    INTEGER, INTENT(IN ) :: nffcst
    INTEGER, INTENT(IN ) :: ifday
    REAL(KIND=r8),    INTENT(IN ) :: tod
    INTEGER, INTENT(IN ) :: idate(4)
    INTEGER, INTENT(IN ) :: idatec(4)
    REAL(KIND=r8),    INTENT(IN ) :: longitude
    REAL(KIND=r8),    INTENT(IN ) :: latitude   
    REAL(KIND=r8),    INTENT(IN ) :: pmand(:)  
    INTEGER, INTENT(IN ) :: it
    CHARACTER(LEN=200), INTENT(IN   ) :: roperm
    CHARACTER(LEN=  7), INTENT(IN   ) :: namef
    CHARACTER(LEN= 10), INTENT(IN   ) :: labeli
    CHARACTER(LEN= 10), INTENT(IN   ) :: labelf
    CHARACTER(LEN=  5), INTENT(IN   ) :: extw
    CHARACTER(LEN=  5), INTENT(IN   ) :: exdw         
    CHARACTER(LEN=  4), INTENT(IN   ) :: trunc
    CHARACTER(LEN=  3), INTENT(IN   ) :: lev
    
    INTEGER :: iyi
    INTEGER :: imi
    INTEGER :: idi
    INTEGER :: ihi
    INTEGER :: iyc
    INTEGER :: imc
    INTEGER :: idc
    INTEGER :: ihc 
    INTEGER :: k
    INTEGER :: ii
    INTEGER :: inv
    INTEGER :: m
    INTEGER :: nn
    INTEGER :: ix
    LOGICAL(KIND=i8) :: inic

    INTEGER,            SAVE :: icall=1
    INTEGER,            SAVE :: is
    CHARACTER(LEN= 10)       :: labelc
    CHARACTER(LEN=  3), SAVE :: ext
    CHARACTER(LEN=  6)       :: extn
    CHARACTER(LEN=  6)       :: exdn
    CHARACTER(LEN= 23)       :: modout
    CHARACTER(LEN= 10)       :: label
    CHARACTER(LEN= 3)        :: cmth(12)
    cmth(1:12)=(/'JAN','FEB','MAR','APR','MAY','JUN',&
                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    inic=(ifday.EQ.0 .AND. tod.EQ.0.0_r8)
    
    IF (icall.EQ.1 .AND. inic) THEN
       ext='icn'
    ELSEIF (icall.EQ.2 .AND. inic) THEN
       icall=3
       ext='inz'
    ELSE
       ext='fct'
    ENDIF
    !modout='/model/dataout/'//trunc//lev//'/'
    modout='/'
    IF (icall .EQ. 1) THEN
       icall=2
       is=INDEX(roperm//' ',' ')-1
       IF (is .LE. 0) is=1
       OPEN(nffcst,file=roperm(1:is)//TRIM(modout)//namef// &
         labeli//extw(1:2)//'fct'//'.'//trunc//lev//'.ctl', &
         status='unknown')
       WRITE(nffcst,'(3A)')'DSET ^'//TRIM(fname(1:17)),'%y4%m2%d2%h2',TRIM(fname(28:40))
       WRITE(nffcst,'(A)')'*'
       WRITE(nffcst,'(A)')'OPTIONS TEMPLATE SEQUENTIAL YREV BIG_ENDIAN'
       WRITE(nffcst,'(A)')'*'
       WRITE(nffcst,'(A)')'UNDEF -2.56E+33'
       WRITE(nffcst,'(A)')'*'
       WRITE(nffcst,'(A)')'TITLE SIMULACAO 1D'
       WRITE(nffcst,'(A)')'*'
       WRITE(nffcst,'(A,I3,A,F8.3,F15.10)') &
                    'XDEF ',1,' LINEAR ',longitude, 1.0_r8
       WRITE(nffcst,'(A,I3,A,F8.3,F15.10)') &
                    'YDEF ',1,' LINEAR ',latitude, 1.0_r8
       WRITE(nffcst,'(A,I3,A)')'ZDEF ',kMaxNew,' LEVELS '
       WRITE(nffcst,'((8X,10F9.6))')(pmand(k),k=1,kMaxNew)
       iyi=idate(4)
       imi=idate(2)
       idi=idate(3)
       ihi=idate(1)
       WRITE(nffcst,'(A,I6,A,I2.2,A,I2.2,A,I4,A)') &
         'TDEF ',it,' LINEAR ',ihi,'Z',idi,cmth(imi),iyi,' 1HR'
       WRITE(nffcst,'(A)')'*'
       WRITE(nffcst,'(A,I3)')'VARS ',mxrq+17
       WRITE(nffcst,'(A)')'TOPO   0 99 '// &
        'TOPOGRAPHY                              (M               )'
       WRITE(nffcst,'(A)')'LSMK   0 99 '// &
        'LAND SEA MASK                           (NO DIM          )'
       WRITE(nffcst,'(A)')'PSLC   0 99 '// &
        'SURFACE PRESSURE                        (Mb              )'
       WRITE(nffcst,'(A,i3,A)')'DIVG ',kMaxNew,' 99 '// &
        'DIVERGENCE                              (1/Sec           )'
       WRITE(nffcst,'(A,i3,A)')'VORT ',kMaxNew,' 99 '// &
        'VORTICITY                               (1/Sec           )'
       WRITE(nffcst,'(A,i3,A)')'UMES ',kMaxNew,' 99 '// &
        'SPECIFIC HUMIDITY                       (No Dim          )'
       WRITE(nffcst,'(A,i3,A)')'TEMP ',kMaxNew,' 99 '// &
        'ABSOLUTE TEMPERATURE                    (K               )'
       WRITE(nffcst,'(A,i3,A)')'UVEL ',kMaxNew,' 99 '// &
        'ZONAL WIND                              (m/Sec           )'
       WRITE(nffcst,'(A,i3,A)')'VVEL ',kMaxNew,' 99 '// &
        'MERIDIONAL WIND                         (m/Sec           )'
       WRITE(nffcst,'(A,i3,A)')'ZORL ', 0  ,' 99 '// &
        'ROUGHNESS LENGTH                        (M               )'
       WRITE(nffcst,'(A,i3,A)')'TSFC ', 0  ,' 99 '// &
        'SURFACE TEMPERATURE                     (K               )'
       WRITE(nffcst,'(A,i3,A)')'TD0S ', 0  ,' 99 '// &
        'deep soil temperature                   (K               )'
       WRITE(nffcst,'(A,i3,A)')'AMDL ', 0  ,' 99 '// &
        'STORAGE ON CANOPY                       (M               )'
       WRITE(nffcst,'(A,i3,A)')'AMSL ', 0  ,' 99 '// &
        'STORAGE ON GROUND                       (M               )'
       WRITE(nffcst,'(A,i3,A)')'USSL ', 0  ,' 99 '// &
        'SOIL WETNESS OF SURFACE                 (No Dim          )'
       WRITE(nffcst,'(A,i3,A)')'UZRS ', 0  ,' 99 '// &
        'SOIL WETNESS OF ROOT ZONE               (No Dim          )'
       WRITE(nffcst,'(A,i3,A)')'UZDS ', 0  ,' 99 '// &
        'SOIL WETNESS OF DRAINAGE ZONE           (No Dim          )'
!*******************************8
    DO m=1,mxavl
      IF (dodia(m) .and. (iavrq(m) > 0)) THEN
        nn=iavrq(m)
        inv=lvrq(nn)
        IF (inv == 1) inv=0
        !WRITE(n,160)reqdg(nn),diag,isg(itavl(m)),lvrq(nn),nurq(nn)
        WRITE(nffcst,'(A,1X,2I3,1X,4A)')alias(nn),inv,99,reqdg(nn),'(',aunits(nurq(nn)),')'        
      END IF
    END DO

    IF(icf.ne.0)THEN
      DO ix=1,icf        
        inv=lvrq(ix)
        IF (inv == 1) inv=0
        !WRITE(n,160)combf(ix),diag,isg(itcf(ix)),lvcf(ix),nucf(ix)
        WRITE(nffcst,'(A,1X,2I3,1X,4A)')alias(ix),inv,99,reqdg(ix),'(',aunits(nurq(ix)),')'        
      END DO
    END IF
!********************************
!       DO ii=1,nof
!        inv=lvrq(ii)
!        IF (inv == 1) inv=0
!        WRITE(nffcst,'(A,1X,2I3,1X,4A)')alias(ii),inv,99,reqdg(ii),'(',aunits(nurq(ii)),')'
!       ENDDO
       WRITE(nffcst,'(A)')'ENDVARS'
       CLOSE(nffcst, STATUS='KEEP')
    ENDIF
    iyi=idate(4)
    imi=idate(2)
    idi=idate(3)
    ihi=idate(1)
    WRITE(label,'(i4.4,3i2.2)')iyi,imi,idi,ihi
    iyc=idatec(4)
    imc=idatec(2)
    idc=idatec(3)
    ihc=idatec(1)
    WRITE(labelc,'(i4.4,3i2.2)')iyc,imc,idc,ihc
    extn(1:2)=extw(1:2)
    extn(3:5)=ext(1:3)
    extn(6:6)='.'
    exdn(1:2)=exdw(1:2)
    IF (ext .EQ. 'icn') THEN
       exdn(3:5)='dic'
    ELSEIF (ext .EQ. 'inz') THEN
       exdn(3:5)='din'
    ELSE
       exdn(3:5)='dir'
    ENDIF
    exdn(6:6)='.'
    !WRITE(*,'(a,3(2x,a))') ' OPNFCT : ',labeli,label,labelc
    CLOSE(nfdiag)
             
    OPEN(nfdiag,file=roperm(1:is)//TRIM(modout)//namef// &
         labeli//labelc//extn//trunc//lev//'.ctl', &
         ACTION='write', STATUS='replace')
 
       WRITE(nfdiag,'(A)')'DSET ^'//TRIM(fname)
       WRITE(nfdiag,'(A)')'*'
       WRITE(nfdiag,'(A)')'OPTIONS SEQUENTIAL YREV'
       WRITE(nfdiag,'(A)')'*'
       WRITE(nfdiag,'(A)')'UNDEF -2.56E+33'
       WRITE(nfdiag,'(A)')'*'
       WRITE(nfdiag,'(A)')'TITLE SIMULACAO 1D'
       WRITE(nfdiag,'(A)')'*'
       WRITE(nfdiag,'(A,I3,A,F8.3,F15.10)') &
                    'XDEF ',1,' LINEAR ',longitude, 1.0_r8
       WRITE(nfdiag,'(A,I3,A,F8.3,F15.10)') &
                    'YDEF ',1,' LINEAR ',latitude, 1.0_r8
       WRITE(nfdiag,'(A,I3,A)')'ZDEF ',kMaxNew,' LEVELS '
       WRITE(nfdiag,'((8X,10F9.6))')(pmand(k),k=1,kMaxNew)
       iyi=idate(4)
       imi=idate(2)
       idi=idate(3)
       ihi=idate(1)
       WRITE(nfdiag,'(A,I3,A,I2.2,A,I2.2,A,I4,A)') &
         'TDEF ',1,' LINEAR ',ihi,'Z',idi,cmth(imi),iyi,' 1HR'
       WRITE(nfdiag,'(A)')'*'
       WRITE(nfdiag,'(A,I3)')'VARS ',mxrq+17
       WRITE(nfdiag,'(A)')'TOPO   0 99 '// &
        'TOPOGRAPHY                              (M               )'
       WRITE(nfdiag,'(A)')'LSMK   0 99 '// &
        'LAND SEA MASK                           (NO DIM          )'
       WRITE(nfdiag,'(A)')'PSLC   0 99 '// &
        'SURFACE PRESSURE                        (Mb              )'
       WRITE(nfdiag,'(A,i3,A)')'DIVG ',kMaxNew,' 99 '// &
        'DIVERGENCE                              (1/Sec           )'
       WRITE(nfdiag,'(A,i3,A)')'VORT ',kMaxNew,' 99 '// &
        'VORTICITY                               (1/Sec           )'
       WRITE(nfdiag,'(A,i3,A)')'UMES ',kMaxNew,' 99 '// &
        'SPECIFIC HUMIDITY                       (No Dim          )'
       WRITE(nfdiag,'(A,i3,A)')'TEMP ',kMaxNew,' 99 '// &
        'ABSOLUTE TEMPERATURE                    (K               )'
       WRITE(nfdiag,'(A,i3,A)')'UVEL ',kMaxNew,' 99 '// &
        'ZONAL WIND                              (m/Sec           )'
       WRITE(nfdiag,'(A,i3,A)')'VVEL ',kMaxNew,' 99 '// &
        'MERIDIONAL WIND                         (m/Sec           )'
       WRITE(nfdiag,'(A,i3,A)')'ZORL ', 0  ,' 99 '// &
        'ROUGHNESS LENGTH                        (M               )'
       WRITE(nfdiag,'(A,i3,A)')'TSFC ', 0  ,' 99 '// &
        'SURFACE TEMPERATURE                     (K               )'
       WRITE(nfdiag,'(A,i3,A)')'TD0S ', 0  ,' 99 '// &
        'deep soil temperature                   (K               )'
       WRITE(nfdiag,'(A,i3,A)')'AMDL ', 0  ,' 99 '// &
        'STORAGE ON CANOPY                       (M               )'
       WRITE(nfdiag,'(A,i3,A)')'AMSL ', 0  ,' 99 '// &
        'STORAGE ON GROUND                       (M               )'
       WRITE(nfdiag,'(A,i3,A)')'USSL ', 0  ,' 99 '// &
        'SOIL WETNESS OF SURFACE                 (No Dim          )'
       WRITE(nfdiag,'(A,i3,A)')'UZRS ', 0  ,' 99 '// &
        'SOIL WETNESS OF ROOT ZONE               (No Dim          )'
       WRITE(nfdiag,'(A,i3,A)')'UZDS ', 0  ,' 99 '// &
        'SOIL WETNESS OF DRAINAGE ZONE           (No Dim          )'
!*******************************8
    DO m=1,mxavl
      IF (dodia(m) .and. (iavrq(m) > 0)) THEN
        nn=iavrq(m)
        inv=lvrq(nn)
        IF (inv == 1) inv=0
        !WRITE(n,160)reqdg(nn),diag,isg(itavl(m)),lvrq(nn),nurq(nn)
        WRITE(nfdiag,'(A,1X,2I3,1X,4A)')alias(nn),inv,99,reqdg(nn),'(',aunits(nurq(nn)),')'
      END IF
    END DO

    IF(icf.ne.0)THEN
      DO ix=1,icf        
        inv=lvrq(ix)
        IF (inv == 1) inv=0
        !WRITE(n,160)combf(ix),diag,isg(itcf(ix)),lvcf(ix),nucf(ix)
        WRITE(nfdiag,'(A,1X,2I3,1X,4A)')alias(ix),inv,99,reqdg(ix),'(',aunits(nurq(ix)),')'
      END DO
    END IF
!********************************
!       DO ii=1,nof
!        inv=lvrq(ii)
!        IF (inv == 1) inv=0
!        WRITE(nfdiag,'(A,1X,2I3,1X,4A)')alias(ii),inv,99,reqdg(ii),'(',aunits(nurq(ii)),')'
!       ENDDO
       WRITE(nfdiag,'(A)')'ENDVARS'
  
  END SUBROUTINE wrfct



  SUBROUTINE opnfct(nfdrct, nfdiag, nffcst, ifday, tod, idate, idatec,&
                    roperm,namef,labeli,labelf,extw,exdw,trunc,lev)
    INTEGER, INTENT(IN ) :: nfdrct
    INTEGER, INTENT(IN ) :: nfdiag
    INTEGER, INTENT(IN ) :: nffcst
    INTEGER, INTENT(IN ) :: ifday
    REAL(KIND=r8),    INTENT(IN ) :: tod
    INTEGER, INTENT(IN ) :: idate(4)
    INTEGER, INTENT(IN ) :: idatec(4)
    
    CHARACTER(LEN=200), INTENT(IN   ) :: roperm
    CHARACTER(LEN=  7), INTENT(IN   ) :: namef
    CHARACTER(LEN= 10), INTENT(IN   ) :: labeli
    CHARACTER(LEN= 10), INTENT(IN   ) :: labelf
    CHARACTER(LEN=  5), INTENT(IN   ) :: extw
    CHARACTER(LEN=  5), INTENT(IN   ) :: exdw         
    CHARACTER(LEN=  4), INTENT(IN   ) :: trunc
    CHARACTER(LEN=  3), INTENT(IN   ) :: lev
    
    INTEGER :: iyi
    INTEGER :: imi
    INTEGER :: idi
    INTEGER :: ihi
    INTEGER :: iyc
    INTEGER :: imc
    INTEGER :: idc
    INTEGER :: ihc
    LOGICAL(KIND=i8) :: inic

    INTEGER,            SAVE :: icall=1
    INTEGER,            SAVE :: is
    CHARACTER(LEN= 10)       :: labelc
    CHARACTER(LEN=  3), SAVE :: ext
    CHARACTER(LEN=  6)       :: extn
    CHARACTER(LEN=  6)       :: exdn
    CHARACTER(LEN= 23)       :: modout
    CHARACTER(LEN= 10)       :: label

    inic=(ifday.EQ.0 .AND. tod.EQ.0.0_r8)
    
    IF (icall.EQ.1 .AND. inic) THEN
       ext='icn'
    ELSEIF (icall.EQ.2 .AND. inic) THEN
       icall=3
       ext='inz'
    ELSE
       ext='fct'
    ENDIF
    !modout='/model/dataout/'//trunc//lev//'/'
    modout='/'
    IF (icall .EQ. 1) THEN
       icall=2
       is=INDEX(roperm//' ',' ')-1
       IF (is .LE. 0) is=1
       OPEN(nffcst,file=roperm(1:is)//TRIM(modout)//namef// &
         labeli//labelf//extw(1:2)//'dir'//'.'//trunc//lev//'.files', &
         status='unknown')
    ENDIF
    iyi=idate(4)
    imi=idate(2)
    idi=idate(3)
    ihi=idate(1)
    WRITE(label,'(i4.4,3i2.2)')iyi,imi,idi,ihi
    iyc=idatec(4)
    imc=idatec(2)
    idc=idatec(3)
    ihc=idatec(1)
    WRITE(labelc,'(i4.4,3i2.2)')iyc,imc,idc,ihc
    extn(1:2)=extw(1:2)
    extn(3:5)=ext(1:3)
    extn(6:6)='.'
    exdn(1:2)=exdw(1:2)
    IF (ext .EQ. 'icn') THEN
       exdn(3:5)='dic'
    ELSEIF (ext .EQ. 'inz') THEN
       exdn(3:5)='din'
    ELSE
       exdn(3:5)='dir'
    ENDIF
    exdn(6:6)='.'
    WRITE(*,'(a,3(2x,a))') ' OPNFCT : ',labeli,label,labelc
    CLOSE(nfdrct)
    CLOSE(nfdiag)
    
    OPEN(nfdrct,file=roperm(1:is)//TRIM(modout)//namef// &
         labeli//labelc//exdn//trunc//lev,form='formatted', &
         ACTION='write', STATUS='replace', IOSTAT=ierr)
    fname=namef//labeli//labelc//extn//trunc//lev
    OPEN(nfdiag,file=roperm(1:is)//TRIM(modout)//namef// &
         labeli//labelc//extn//trunc//lev,form='unformatted', &
         ACTION='write', STATUS='replace', IOSTAT=ierr)
    WRITE(nffcst,'(a)')roperm(1:is)//TRIM(modout)//namef// &
         labeli//labelc//exdn//trunc//lev
    WRITE(nffcst,'(a)')roperm(1:is)//TRIM(modout)//namef// &
         labeli//labelc//extn//trunc//lev         
  
  END SUBROUTINE opnfct







  SUBROUTINE wrprog (iudrct,iudiag,ifday ,tod   ,idate ,idatec,qrot  , &
                     qdiv  ,qq    ,qlnp  ,qtmp  ,zorl  ,gtsea ,td0   , &
                     capac0,w0    ,imask ,nexp  ,jttl  ,iufcst,del   , &
                     qgzs  ,lsmk  ,ijMaxGauQua  ,kmax  ,imax  ,jmax  , &
                     roperm,namef ,labeli,labelf,extw  ,exdw  ,trunc , &
                     lev,longitude,latitude,pmand,it,gu  ,gv   )
    INTEGER           , INTENT(IN   ) :: iudrct
    INTEGER           , INTENT(IN   ) :: iudiag
    INTEGER           , INTENT(IN   ) :: ifday
    INTEGER           , INTENT(IN   ) :: ijMaxGauQua
    INTEGER           , INTENT(IN   ) :: imax
    INTEGER           , INTENT(IN   ) :: jmax
    INTEGER           , INTENT(IN   ) :: kmax
    REAL(KIND=r8)              , INTENT(IN   ) :: tod
    INTEGER           , INTENT(IN   ) :: idate (:)
    INTEGER           , INTENT(IN   ) :: idatec(:) 
    REAL(KIND=r8)              , INTENT(IN   ) :: qlnp  (iMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: qtmp  (iMax,kMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: qdiv  (iMax,kMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: qrot  (iMax,kMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: qq    (iMax,kMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: gu    (iMax,kMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: gv    (iMax,kMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: lsmk  (iMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: zorl  (iMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: gtsea (iMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: td0   (iMax,jMax)
    REAL(KIND=r8)              , INTENT(IN   ) :: capac0(:,:,:)
    REAL(KIND=r8)              , INTENT(IN   ) :: w0    (:,:,:)
    INTEGER(KIND=i8)           , INTENT(IN   ) :: imask (iMax,jMax)
    CHARACTER(LEN= 4) , INTENT(IN   ) :: nexp
    CHARACTER(LEN=40) , INTENT(IN   ) :: jttl
    INTEGER           , INTENT(IN   ) :: iufcst
    REAL(KIND=r8)              , INTENT(IN   ) :: del   (kMaxNew)
    REAL(KIND=r8)              , INTENT(IN   ) :: qgzs  (:,:)
    CHARACTER(LEN=200), INTENT(IN   ) :: roperm
    CHARACTER(LEN=  7), INTENT(IN   ) :: namef
    CHARACTER(LEN= 10), INTENT(IN   ) :: labeli
    CHARACTER(LEN= 10), INTENT(IN   ) :: labelf
    CHARACTER(LEN=  5), INTENT(IN   ) :: extw
    CHARACTER(LEN=  5), INTENT(IN   ) :: exdw         
    CHARACTER(LEN=  4), INTENT(IN   ) :: trunc
    CHARACTER(LEN=  3), INTENT(IN   ) :: lev
    REAL(KIND=r8)              , INTENT(IN   ) :: longitude
    REAL(KIND=r8)              , INTENT(IN   ) :: latitude    
    REAL(KIND=r8)              , INTENT(IN   ) :: pmand(:)
    INTEGER           , INTENT(IN   ) :: it
    REAL(KIND=r8)                              :: work(iMax,jMax)
    REAL(KIND=r8)                              :: work_in(iMax,jMax)
    INTEGER                           :: k,ij,j,i
    INTEGER                           :: nlev
    INTEGER                           :: ihdim
    REAL(KIND=r8)                              :: confac
    INTEGER                           :: ifday4
    INTEGER                           :: idat4(4)
    INTEGER                           :: idat4c(4)
    REAL(KIND=r8)                              :: tod4
    INTEGER                           :: outctl=150
        
    LOGICAL(KIND=i8)               , PARAMETER :: toCol=.FALSE.
    REAL(KIND=r8)                              :: aux1(ijMaxGauQua)
    REAL(KIND=r8)                              :: aux2(ijMaxGauQua,kMax)
    REAL(KIND=r8)                              :: aux3(ijMaxGauQua)
    
    CALL opnfct(iudrct,iudiag,iufcst,ifday,tod,idate,idatec,&
                roperm,namef,labeli,labelf,extw,exdw,trunc,lev)
    CALL wrfct(iudrct,71,72,ifday,tod,idate,idatec,&
               roperm,namef,labeli,labelf,extw,exdw,trunc,lev,longitude,latitude,&
               pmand,it)
    CALL WriteDir(iudrct, idate,idatec(1), idatec(3),idatec(2),&
                  idatec(4),del,tod)    
    ifday4=ifday
    tod4=tod
    DO k=1,4
       idat4(k)=idate(k)
       idat4c(k)=idatec(k)
    ENDDO
    
    !CALL WriteProgHead(iudiag, ifday4, tod4, idat4, idat4c)
    !
    !     topography
    !
    !     
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux1(ij)=qgzs(i,j) 
          ij=ij+1
      END DO
    END DO   
    CALL WriteField(iudiag, aux1)
    !
    !     land sea mask
    !
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux1(ij)=lsmk(i,j) 
          ij=ij+1
      END DO
    END DO    
    CALL WriteField(iudiag, aux1)
    !          
    !          ln surface pressure
    !     
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux1(ij)=exp(qlnp(i,j) )
          ij=ij+1
      END DO
    END DO   
    CALL WriteField(iudiag, aux1)
    !
    !     divergence 
    !     
    DO k=1,kMax
      ij=1
      DO j=1,jMax
         DO i=1,iMax
           aux2(ij,k)=qdiv(i,k,j) 
           ij=ij+1
         END DO
      END DO   
    END DO
    CALL WriteField(iudiag, aux2)
    !
    !     vorticity
    !
    DO k=1,kMax
      ij=1
      DO j=1,jMax
         DO i=1,iMax
           aux2(ij,k)=qrot(i,k,j) 
           ij=ij+1
         END DO
      END DO   
    END DO
    CALL WriteField(iudiag, aux2)    
    !
    !     specific humidity
    !     
    DO k=1,kMax
      ij=1
      DO j=1,jMax
         DO i=1,iMax
           aux2(ij,k)=qq(i,k,j) 
           ij=ij+1
         END DO
      END DO   
    END DO
    CALL WriteField(iudiag, aux2)
    !     
    !     virtual temperature
    !     
    DO k=1,kMax
      ij=1
      DO j=1,jMax
         DO i=1,iMax
           aux2(ij,k)=qtmp(i,k,j) + tov(k)
           ij=ij+1
         END DO
      END DO   
    END DO
    CALL WriteField(iudiag, aux2)
    !     
    !    Componente da Velocidade Zonal
    !     
    DO k=1,kMax
      ij=1
      DO j=1,jMax
         DO i=1,iMax
           aux2(ij,k)=gu(i,k,j)
           ij=ij+1
         END DO
      END DO   
    END DO
    CALL WriteField(iudiag, aux2)
    !     
    !   Componente da Velocidade Meridional
    !     
    DO k=1,kMax
      ij=1
      DO j=1,jMax
         DO i=1,iMax
           aux2(ij,k)=gv(i,k,j)
           ij=ij+1
         END DO
      END DO   
    END DO
    CALL WriteField(iudiag, aux2)
    !     
    !     surface roughness 
    !     
    ij=1
    work_in=zorl
    DO j=1,jMax
       DO i=1,iMax
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    !
    !     surface temperature
    !     
    ij=1
    work_in=gtsea
    DO j=1,jMax
       DO i=1,iMax
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    !
    !     deep soil temperature 
    !     
    CALL reord (td0,    1, work, 1, imask, gtsea, 'TD  ')
    work_in=work
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              PRINT*,'td0=',td0
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    !     
    !     storage on canopy
    !     
    CALL reord (capac0, 2, work, 1, imask, gtsea, 'CAPC') 
    work_in=work
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    !     
    !     storage on ground cover
    !     
    CALL reord (capac0, 2, work, 2, imask, gtsea, 'CAPG')
    work_in=work
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    !
    !     wetness of surface zone
    !     
    CALL reord (w0,     3, work, 1, imask, gtsea, 'W1  ')
    work_in=work
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    !     
    !     wetness of root zone
    !     
    CALL reord (w0,     3, work, 2, imask, gtsea, 'W2  ')
    work_in=work
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    !     
    !     wetness of drainage zone
    !     
    CALL reord (w0,     3, work, 3, imask, gtsea, 'W3  ')
    work_in=work
    ij=1
    DO j=1,jMax
       DO i=1,iMax
              aux3(ij)=work_in(i,j) 
          ij=ij+1
      END DO
    END DO
    CALL WriteField(iudiag, aux3)
    
    IF(ifprt(95).GE.1)WRITE(nfprt,5000)idate,ifday,tod,idatec
5000 FORMAT(' DONE WITH WRPROG. MODEL STARTED ',3I3,I5/' NOW AT',I8, &
         ' DAYS AND',F8.1,' SECONDS.  CURRENT DATE IS',3I3,I5)
  END SUBROUTINE wrprog
END MODULE Diagnostics
