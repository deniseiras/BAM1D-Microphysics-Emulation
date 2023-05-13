MODULE Sfc_Ibis_Interface
  !
  ! InitFieldsIbis__vegin
  !              |                
  !              |__RD_PARAM
  !              |
  !              |__readit__ibismap__cellbox
  !              |                |
  !              |                |__cellbox 
  !              |
  !              |__climanl__existence
  !              |
  !              |__initial__coldstart
  !              |        |
  !              |        |__restart__existence
  !              |        |
  !              |        |__inisurf
  !              |        |
  !              |        |__inisnow
  !              |        |
  !              |        |__inisoil
  !              |        |    
  !              |        |__iniveg
  !              |        | 
  !              |        |__inisum
  !
  ! IbisDrv__Ibis__pheno
  !       |     |
  !       |     |__soilbgc
  !       |     |
  !       |     |__lsxmainn_______setsoi
  !       |     |              |
  !       |     |              |__fwetcal
  !       |     |              |
  !       |     |              |__solset
  !       |     |              |
  !       |     |              |__solsur
  !       |     |              |
  !       |     |              |__solalb___twostr__twoset
  !       |     |              |        |
  !       |     |              |        |__twostr__twoset
  !       |     |              |
  !       |     |              |__solarf
  !       |     |              |
  !       |     |              |__irrad
  !       |     |              |
  !       |     |              |__cascade___mix
  !       |     |              |         | 
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__steph2o__mix
  !       |     |              |         |         |
  !       |     |              |         |         |__mix
  !       |     |              |         |
  !       |     |              |         |__steph2o__mix
  !       |     |              |         |         |
  !       |     |              |         |         |__mix
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__steph2o__mix
  !       |     |              |         |         |
  !       |     |              |         |         |__mix  
  !       |     |              |         |
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              |         |
  !       |     |              |         |__mix
  !       |     |              | 
  !       |     |              |__fwetcal
  !       |     |              |
  !       |     |              |__canopy__canini
  !       |     |              |       |
  !       |     |              |       |__drystress
  !       |     |              |       |
  !       |     |              |       |__turcof__fstrat
  !       |     |              |       |       |
  !       |     |              |       |       |__fstrat
  !       |     |              |       |
  !       |     |              |       |__stomata
  !       |     |              |       |
  !       |     |              |       |__turvap___impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp
  !       |     |              |       |       |
  !       |     |              |       |       |__impexp2
  !       |     |              |       |       |
  !       |     |              |       |       |__linsolve
  !       |     |              |       |
  !       |     |              |       |__tscreen
  !       |     |              |       |
  !       |     |              |       |__tscreen
  !       |     |              |
  !       |     |              |__cascad2__steph2o2
  !       |     |              |        |
  !       |     |              |        |__steph2o2
  !       |     |              |        |
  !       |     |              |        |__steph2o2
  !       |     |              |
  !       |     |              |__noveg
  !       |     |              |
  !       |     |              |__snow__snowheat___tridia
  !       |     |              |     |
  !       |     |              |     |__vadapt
  !       |     |              |     |
  !       |     |              |     |__MixSnow
  !       |     |              |
  !       |     |              |__soilctl__soilh2o__tridia
  !       |     |                       |
  !       |     |                       |__soilheat__tridia
  !       |     |                       |
  !       |     |                       |__wadjust
  !       |     |                       |
  !       |     |                       |__wadjust
  !       |     |
  !       |     |
  !       |     |__sumnow
  !       |     |
  !       |     |__sumday
  !       |     |
  !       |     |__summonth
  !       |     |
  !       |     |__sumyear
  !       |     |
  !       |     |__solset
  !       |     |
  !       |     |__solsur
  !       |     |
  !       |     |__solalb___twostr__twoset
  !       |     |       |
  !       |     |       |___twostr__twoset
  !       |     |
  !       |     |__solarf
  !       |
  !       |__co2
  !       |
  !       |__dynaveg2__fire
  !       |         |
  !       |         |__vegmap
  !       |
  !       |__climanl2__existence
  !
  !
  ! Albedo_IBIS__fwetcal
  !           |
  !           |__solset
  !           |
  !           |__solsur
  !           |
  !           |__solalb___twostr__twoset
  !                   |
  !                   |___twostr__twoset

  USE Constants, ONLY :     &
       oceald   ,     &
       icealn   ,     &
       icealv   ,     &
       r8,i8,r4,i4

  USE Sfc_Ibis_Fiels , ONLY : spinmax ,isimco2 ,doalb,  &
       isimveg ,isimfire,nband, nsoilay, nsnolay,nlpoints,&
       pi,stef,vonk,grav,tmelt,hfus,hvap,hsub,ch2o,npft, &
       cice,cair,cvap,rair,rvap,cappa,rhow,co2init,epsilon,&
       ginvap  ,gsuvap  ,gtrans,gtransu,gtransl,grunof,& 
       gdrain  ,gadjust ,a10td,a10ancub,a10ancuc,a10ancls, & 
       a10ancl3,a10ancl4,a10scalparamu,a10daylightu,a10scalparaml,&
       a10daylightl,vmax_pft,tau15,kc15,ko15,cimax,gammaub,alpha3,&
       theta3,beta3,coefmub,coefbub,gsubmin,gammauc,coefmuc,coefbuc ,&
       gsucmin,gammals,coefmls,coefbls,gslsmin,gammal3,coefml3, &
       coefbl3,gsl3min,gammal4,alpha4,theta4,beta4,coefml4,coefbl4, &
       gsl4min, wliqu,wliqumax,wsnou,wsnoumax,tum,tu,tu0,wliqs,wliqsmax,& 
       wsnos,wsnosmax,ts,wliql,wliqlmax,wsnol,wsnolmax,tl,topparu,&
       topparl,fl,fu,lai,sai,rhoveg,tauveg,orieh,oriev,wliqmin, &
       wsnomin,t12,tdripu,tblowu,tdrips,tblows,t34,tdripl,tblowl,&
       ztop,alaiml,zbot,alaimu,froot,q34,q12,su,cleaf,dleaf,ss        ,& 
       cstem,dstem,sl,cgrass,ciub,ciuc,exist,csub,gsub,csuc,gsuc,&
       agcub,agcuc,ancub,ancuc,totcondub,totconduc,cils,cil3,&         
       cil4,csls,gsls,csl3,gsl3,csl4,gsl4,agcls,agcl4,agcl3,ancls,&         
       ancl4,ancl3,totcondls,totcondl3,totcondl4,chu,chs,chl,frac,& 
       tlsub,z0sno,rhos,consno,hsnotop,hsnomin,fimin,fimax,fi,&     
       tsno,hsno,sand,clay,poros,wsoim,wsoi,wsoi0,wisoi,consoi,zwpmax, wpud,&        
       wipud,wpudmax, qglif ,tsoim ,tsoi,tsoi0,hvasug,hvasui,albsav,albsan,&
       tg,ti,z0soi,swilt,sfield,stressl,stressu,stresstl,stresstu,&
       csoi,rhosoi,hsoi,suction,bex,upsoiu,upsoil,heatg,heati,&
       hydraul,porosflo,ibex,bperm,hflo,o2conc,co2conc,cbiow,&
       sapfrac,cbior ,adrain,adsnow,adaet,adtrunoff,& 
       adsrunoff,addrainage,adrh,adsnod,adsnof,adwsoi,adtsoi, &
       adwisoi,adtlaysoi,adwlaysoi,adwsoic,adtsoic,adco2mic,adco2root, & 
       adco2soi,adco2ratio,adnmintot,decompl,decomps,tnmin,amrain,&     
       amsnow,amaet,amtrunoff,amsrunoff,amdrainage,amtemp,&    
       amqa,amsolar,amirup,amirdown,amsens,amlatent,amlaiu,& 
       amlail,amtsoi,amwsoi,amwisoi,amvwc,amawc,amsnod,amsnof,amnpp,& 
       amnpptot,amco2mic,amco2root,amco2soi,amco2ratio,amneetot, &
       amnmintot,amts2,amtransu,amtransl,amsuvap,aminvap,amalbedo,&
       amtsoil,amwsoil,amwisoil,tnpptot,&
       aysolar,ayirup,ayirdown,aysens,aylatent,ayprcp,&    
       ayaet,aytrans,aytrunoff,aysrunoff,aydrainage,aydwtot,aywsoi, &
       aywisoi,aytsoi,ayvwc,ayawc,aystresstu,aystresstl,aygpp,aygpptot,&  
       aynpp,aynpptot,ayco2mic,ayco2root,ayco2soi,ayneetot,ayrootbio, & 
       aynmintot,ayalit,ayblit,aycsoi,aycmic,ayanlit,aybnlit,aynsoi,ayalbedo, &    
       totalit,totrlit,totcsoi,totcmic,totanlit,totrnlit,totnsoi,&   
       totnmic,totlit,totfall,totnlit,firefac,wtot,storedn,yrleach,& 
       ynleach,falll,fallr,fallw,clitlm,clitls,clitrm,clitrs,clitwm,&  
       clitws,csoislop,csoislon,csoipas,clitll,clitrl,clitwl,tw,&
       tc,agddu,tempu,agddl,templ,dropu,dropls,dropl4,dropl3,plai,&
       deltat,gdd0,gdd0this,tcthis,twthis,tcmin,gdd5,gdd5this,TminL,&        
       TminU,Twarm,GDD,aleaf,awood,cbiol,aroot,specla,td ,vzero,&  
       biomass,totlaiu,totlail,totbiou,totbiol,woodnorm,vegtype0,&
       tauwood0,tauwood,tauleaf,tauroot,xminlai,cdisturb,ayanpp,&
       ayanpptot,asurd,asuri,ndaypm,idateprev,iyear0,&
       ndtimes,nmtimes,nytimes,nppdummy,tco2root,&
       tneetot,tco2mic,iMaskIBIS

  USE Sfc_Ibis_LsxMain  , ONLY : co2,lsxmain, fwetcal,solset,solsur ,solalb,solarf

  USE Sfc_Ibis_Vegetation, ONLY : dynaveg1,dynaveg2,climanl2,gdiag,vdiag,&
       pheno ,sumnow,sumday,summonth,sumyear,soilbgc

  USE FieldsPhysics, ONLY: Tsfc0,Qsfc0,Tsfcm,Qsfcm,w0,wm,capac0,capacm,td0,tdm,tcm,tc0,tgm,tg0

  USE module_sf_ocean, ONLY: sf_gfs_seaice_wrapper,OCEANML

  USE OceanModel         , Only : SF_EXCH

  USE Parallelism, ONLY: &
       MsgOne,           &
       myId,             &
       maxNodes
  USE Diagnostics, ONLY: &
       updia, &
       dodia, &
       nDiag_biomau, &
       nDiag_biomal, &
       nDiag_tlaiup, & 
       nDiag_tlailw, & 
       nDiag_tstnsp, &
       nDiag_wsttot, &
       nDiag_lidecf, &
       nDiag_somdfa, &
       nDiag_facuca, &
       nDiag_fsfclc, &
       nDiag_frsnow, &
       nDiag_insnpp, & ! instantaneous npp (mol-CO2 / m-2 / second)
       nDiag_insnee, & ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
       nDiag_grbdy0, & ! annual total growing degree days for current year > 0C
       nDiag_grbdy5, & ! annual total growing degree days for current year > 5C  
       nDiag_avet2m, & ! monthly average 2-m surface-air temperature 
       nDiag_monnpp, & ! monthly total npp for ecosystem (kg-C/m**2/month)
       nDiag_monnee, & ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
       nDiag_yeanpp, & ! annual total npp for ecosystem (kg-c/m**2/yr)
       nDiag_yeanee, & ! annual total NEE for ecosystem (kg-C/m**2/yr) 
       nDiag_upclai, & ! upper canopy single-sided leaf area index (area leaf/area veg)
       nDiag_lwclai, & ! lower canopy single-sided leaf area index (area leaf/area veg)
       nDiag_pfts01, & ! pft tropical broadleaf evergreen trees
       nDiag_pfts02, & ! pft tropical broadleaf drought-deciduous trees
       nDiag_pfts03, & ! pft warm-temperate broadleaf evergreen trees
       nDiag_pfts04, & ! pft temperate conifer evergreen trees
       nDiag_pfts05, & ! pft temperate broadleaf cold-deciduous trees
       nDiag_pfts06, & ! pft boreal conifer evergreen trees
       nDiag_pfts07, & ! pft boreal broadleaf cold-deciduous trees
       nDiag_pfts08, & ! pft boreal conifer cold-deciduous trees
       nDiag_pfts09, & ! pft evergreen shrubs
       nDiag_pfts10, & ! pft cold-deciduous shrubs
       nDiag_pfts11, & ! pft warm (c4) grasses
       nDiag_pfts12, & ! pft cool (c3) grasses
       nDiag_biol01, & ! cbiol tropical broadleaf evergreen trees
       nDiag_biol02, & ! cbiol tropical broadleaf drought-deciduous trees
       nDiag_biol03, & ! cbiol warm-temperate broadleaf evergreen trees
       nDiag_biol04, & ! cbiol temperate conifer evergreen trees
       nDiag_biol05, & ! cbiol temperate broadleaf cold-deciduous trees
       nDiag_biol06, & ! cbiol boreal conifer evergreen trees
       nDiag_biol07, & ! cbiol boreal broadleaf cold-deciduous trees
       nDiag_biol08, & ! cbiol boreal conifer cold-deciduous trees
       nDiag_biol09, & ! cbiol evergreen shrubs
       nDiag_biol10, & ! cbiol cold-deciduous shrubs
       nDiag_biol11, & ! cbiol warm (c4) grasses
       nDiag_biol12, & ! cbiol cool (c3) grasses
       nDiag_ynpp01, & ! ynpp tropical broadleaf evergreen trees
       nDiag_ynpp02, & ! ynpp tropical broadleaf drought-deciduous trees
       nDiag_ynpp03, & ! ynpp warm-temperate broadleaf evergreen trees
       nDiag_ynpp04, & ! ynpp temperate conifer evergreen trees
       nDiag_ynpp05, & ! ynpp temperate broadleaf cold-deciduous trees
       nDiag_ynpp06, & ! ynpp boreal conifer evergreen trees
       nDiag_ynpp07, & ! ynpp boreal broadleaf cold-deciduous trees
       nDiag_ynpp08, & ! ynpp boreal conifer cold-deciduous trees
       nDiag_ynpp09, & ! ynpp evergreen shrubs
       nDiag_ynpp10, & ! ynpp cold-deciduous shrubs
       nDiag_ynpp11, & ! ynpp warm (c4) grasses
       nDiag_ynpp12, & ! ynpp cool (c3) grasses
       nDiag_cmontp, & !coldest monthly temperature                         (C)
       nDiag_wmontp, & !warmest monthly temperature                         (C)
       nDiag_atogpp, & !annual total gpp for ecosystem                             (kg-c/m**2/yr)
       nDiag_toigpp, & !instantaneous gpp                           (mol-CO2 / m-2 / second)
       nDiag_fxcsol, & !instantaneous fine co2 flux from soil                           (mol-CO2 / m-2 / second)
       nDiag_mcsoil, & !instantaneous microbial co2 flux from soil         (mol-CO2 / m-2 / second)
       nDiag_cagcub, & !canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
       nDiag_cagcuc, & !canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
       nDiag_cagcls, & !canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
       nDiag_cagcl4, & !canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
       nDiag_cagcl3, & !canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
       nDiag_cancub, & !canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
       nDiag_cancuc, & !canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
       nDiag_cancls, & !canopy average net photosynthesis rate - shrubs         (mol_co2 m-2 s-1)
       nDiag_cancl4, & !canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
       nDiag_cancl3, & !canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
       nDiag_cicoub, & !intercellular co2 concentration - broadleaf        (mol_co2/mol_air)
       nDiag_cicouc, & !intercellular co2 concentration - conifer        (mol_co2/mol_air)
       nDiag_cscoub, & !leaf boundary layer co2 concentration - broadleaf     (mol_co2/mol_air)
       nDiag_gscoub, & !upper canopy stomatal conductance - broadleaf                (mol_co2 m-2 s-1)
       nDiag_cscouc, & !leaf boundary layer co2 concentration - conifer         (mol_co2/mol_air)
       nDiag_gscouc, & !upper canopy stomatal conductance - conifer        (mol_co2 m-2 s-1)
       nDiag_cicols, & !intercellular co2 concentration - shrubs        (mol_co2/mol_air)
       nDiag_cicol3, & !intercellular co2 concentration - c3 plants        (mol_co2/mol_air)
       nDiag_cicol4, & !intercellular co2 concentration - c4 plants        (mol_co2/mol_air)
       nDiag_cscols, & !leaf boundary layer co2 concentration - shrubs          (mol_co2/mol_air)
       nDiag_gscols, & !lower canopy stomatal conductance - shrubs        (mol_co2 m-2 s-1)
       nDiag_cscol3, & !leaf boundary layer co2 concentration - c3 plants     (mol_co2/mol_air)
       nDiag_gscol3, & !lower canopy stomatal conductance - c3 grasses          (mol_co2 m-2 s-1)
       nDiag_cscol4, & !leaf boundary layer co2 concentration - c4 plants     (mol_co2/mol_air)
       nDiag_gscol4, & !lower canopy stomatal conductance - c4 grasses          (mol_co2 m-2 s-1)
       nDiag_tcthis, & !coldest monthly temperature of current year          (C)
       nDiag_twthis    !warmest monthly temperature of current year          (C)

  USE Options, ONLY: OCFLUX,omlmodel,oml_hml0

  IMPLICIT NONE
  PRIVATE

  REAL (KIND=r8), PARAMETER   :: z0ice  =       0.001e0_r8! 
  REAL(KIND=r8)   , PARAMETER   :: rgas   = 287.0_r8!  dry air gas constant (J deg^-1 kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: kapa   = 0.2861328125_r8!  rgas/cp (unitless)
  REAL(KIND=r8)   , PARAMETER   :: cp     = rgas/kapa!  specific of heat of dry air at constant pressure (J deg^-1 kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: hltm   = 2.52e6_r8            !  latent heat of vaporization (J kg^-1)
  REAL(KIND=r8)   , PARAMETER   :: stefan = 5.67e-8_r8           !  Stefan-Boltzmann constant
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_AIRMW  = 28.97_r8                  ! kg/Kmole
  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_H2OMW  = 18.01_r8                  ! kg/Kmole

  REAL(KIND=r8), PARAMETER, PUBLIC :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0_r8   ! --

  REAL(KIND=r8)   :: rbyg

  PUBLIC :: Ibis_Interface,Albedo_IBIS 

CONTAINS

  SUBROUTINE Ibis_Interface(intg                ,istrt             , &
       jdt                ,jb                   ,dtc3x             ,nCols             , &
       ktm                ,initlz               , &
       kt                 ,iswrad               ,ilwrad            ,kMax              , &
       tod                ,idatec               ,filta             ,epsflt            , &
       gt                 ,gq                   ,gu                ,gv                , &
       gps                ,tmtx                 ,qmtx              ,umtx              , &
       zenith             ,colrad               ,fira2             ,xvisb             , &
       xvisd              ,xnirb                ,xnird             ,ppli              , &
       ppci               ,snow                 ,sigki             ,delsig            , & 
       tseam              ,tsea                 ,mskant            ,speedm            , &
       slrad              ,tsurf                ,qsurf             ,zorl              , &
       taux               ,tauy                 ,sens              ,evap              , & 
       umom               ,vmom                 ,rmi               ,rhi               , &
       z0                 ,ustar                ,hc                ,hg                , & 
       ec                 ,eg                   ,ts2               ,qs2               , &
       qsfc               ,tsfc                 ,z0sea             ,d                 , &
       cu                 ,imask                ,Ustarm            ,tgrd              , &
       roff               ,ect                  ,eci               ,egt               , &
       egi                ,egs                  ,rho               ,bstar             , &
       HML                ,HUML                 ,HVML              ,TSK               , &
       cldtot             ,ySwSfcNet            ,LwSfcNet          ,pblh              , &
       QCF                ,QCL)

    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: jdt
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    INTEGER, INTENT(IN   ) :: idatec (4)  
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: jb
    INTEGER, INTENT(IN   ) :: ktm
    INTEGER, INTENT(IN   ) :: initlz
    INTEGER, INTENT(IN   ) :: kt
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8)   , INTENT(IN   ) :: dtc3x         ! model timestep (seconds)
    REAL(KIND=r8)   , INTENT(IN   ) :: gt     (nCols,kMax)          ! air temperature (K)
    REAL(KIND=r8)   , INTENT(IN   ) :: gq     (nCols,kMax)          ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8)   , INTENT(IN   ) :: gu     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gv     (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: gps    (nCols)          ! surface pressure (nPa)
    REAL(KIND=r8)   , INTENT(IN   ) :: xvisb  (nCols)!solad(i,1).Downward Surface shortwave fluxe visible beam    (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: xvisd  (nCols)!solai(i,1) !.Downward Surface shortwave fluxe visible diffuse (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: xnirb  (nCols)!solad(i,2).Downward Surface shortwave fluxe Near-IR beam    (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: xnird  (nCols)!solai(1,2) !.Downward Surface shortwave fluxe Near-IR diffuse (cloudy)
    REAL(KIND=r8)   , INTENT(IN   ) :: fira2  (nCols)          ! incoming ir flux (W m-2)
    REAL(KIND=r8)   , INTENT(IN   ) :: zenith (nCols)          ! cosine of solar zenith angle
    REAL(KIND=r8)   , INTENT(IN   ) :: colrad (nCols)           
    REAL(KIND=r8)   , INTENT(IN   ) :: ppli   (nCols)! Precipitation rate ( large scale )       (mm/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: ppci   (nCols)! Precipitation rate ( cumulus )           (mm/s)
    REAL(KIND=r8)   , INTENT(IN   ) :: snow   (nCols)! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8)   , INTENT(IN   ) :: sigki  
    REAL(KIND=r8)   , INTENT(IN   ) :: delsig   

    REAL(KIND=r8)   , INTENT(INOUT) :: tseam  (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsea   (nCols)

    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(ncols)

    REAL(KIND=r8)   , INTENT(INOUT) :: speedm (nCols) ! wind speed (m s-1)
    REAL(KIND=r8)   , INTENT(INOUT) :: tmtx   (nCols,3)
    REAL(KIND=r8)   , INTENT(INOUT) :: qmtx   (nCols,3)
    REAL(KIND=r8)   , INTENT(INOUT) :: umtx   (nCols,4)
    REAL(KIND=r8)   , INTENT(IN   ) :: slrad  (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsurf  (nCols) 
    REAL(KIND=r8)   , INTENT(IN   ) :: qsurf  (nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: zorl   (nCols) 
    REAL(KIND=r8)   , INTENT(OUT  ) :: taux   (nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: tauy   (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: ts2    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: qs2    (nCols)
    REAL(KIND=r8)   , INTENT(IN   ) :: filta
    REAL(KIND=r8)   , INTENT(IN   ) :: epsflt
    INTEGER         , INTENT(IN   ) :: intg
    INTEGER         , INTENT(IN   ) :: istrt
    CHARACTER(len=*), INTENT(IN   ) :: iswrad
    CHARACTER(len=*), INTENT(IN   ) :: ilwrad
    REAL(KIND=r8)   , INTENT(INOUT) :: sens(1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: evap(1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: umom  (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: vmom  (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: rmi   (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: rhi   (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tsfc  (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: qsfc  (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: z0    (1:nCols) 
    REAL(KIND=r8)   , INTENT(INOUT) :: ustar (1:nCols )
    REAL(KIND=r8)   , INTENT(OUT  ) :: hc (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: hg (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: ec (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: eg (1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: z0sea (1:nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: d(1:nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: cu(1:nCols)
    INTEGER(KIND=i8), INTENT(IN   ) :: imask    (nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: Ustarm(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: tgrd(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: roff(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: ect(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: eci(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: egt(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: egi(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: egs(nCols)
    REAL(KIND=r8)   , INTENT(INOUT) :: rho(nCols)
    REAL(KIND=r8)   , INTENT(OUT  ) :: bstar(nCols)

    REAL(KIND=r8),    INTENT(INOUT) :: HML  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HUML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HVML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: TSK  (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: cldtot (nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: ySwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: LwSfcNet (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: pblh (ncols)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCF(nCols,kMax)
    REAL(KIND=r8)   , INTENT(IN   ) :: QCL(nCols,kMax)
    !  idatec(1)....hour(00/12)
    !  idatec(2)....month
    !  idatec(3)....day of month
    !  idatec(4)....year
    !   INCLUDE 'comatm.h'
    REAL(KIND=r8)   :: cond  (1:nCols)
    REAL(KIND=r8)   :: stor  (1:nCols)
    REAL(KIND=r8)   :: rnet  (1:nCols)
    REAL(KIND=r8)         :: ta     (nCols)        ! air temperature (K)
    REAL(KIND=r8)         :: qa     (nCols)        ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8)         :: ux     (nCols)
    REAL(KIND=r8)         :: uy     (nCols)
    REAL(KIND=r8)         :: psurf  (nCols)        ! surface pressure (Pa)
    REAL(KIND=r8)         :: solad  (nCols,nband) ! direct downward solar flux (W m-2)
    REAL(KIND=r8)         :: solai  (nCols,nband) ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8)         :: fira   (nCols)        ! incoming ir flux (W m-2)
    REAL(KIND=r8)         :: coszen (nCols)        ! cosine of solar zenith angle
    REAL(KIND=r8)         :: raina  (nCols)        ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8)         :: snowa  (nCols)        ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8)         :: ua     (nCols)        ! wind speed (m s-1)
    REAL(KIND=r8)         :: tmin   (nCols)        ! minimun soil temperature (K)
    REAL(KIND=r8)         :: tmax   (nCols)        ! maximun soil temperature (K)
    INTEGER, PARAMETER :: ndaypy=365  ! number of days per year(365or366)
    REAL(KIND=r8)         :: tgpptot(nCols) 
    REAL(KIND=r8)         :: disturbf(nCols)    ! local annual fire disturbance regime (m2/m2/yr)
    REAL(KIND=r8)         :: disturbo(nCols)    ! local fraction of biomass pool lost every year to disturbances other than fire
    REAL(KIND=r8)   :: fvapa    (nCols)         ! local ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8)   :: fsena    (nCols)         ! local ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8)   :: xsea       (1:nCols) 
    REAL(KIND=r8)   :: bstar1       (1:nCols)
    REAL(KIND=r8)   :: diag       (1:nCols)
    REAL(KIND=r8)   :: firb     (nCols)          ! local ! net upward ir radiation at reference
    REAL(KIND=r8)   :: GSW (nCols)
    REAL(KIND=r8)   :: GLW (nCols)

    ! atmospheric level za (W m-2)

    !
    ! Arguments (input)
    !
    LOGICAL :: InitMod
    !REAL(KIND=r8)    :: calday     ! current julian day (1-365.99)
    !  INTEGER :: iday     ! current day in month (1-31 or 1-30, or 1-28)
    INTEGER :: iday         ! day number  (passed in)
    !INTEGER :: idayprev   ! day in month of previous timestep

    !   integer:: imonth         ! current month (1 - 12)
    INTEGER :: imonth         ! month number (passed in)   
    INTEGER :: imonthprev ! month of previous timestep 

    INTEGER :: iyear
    !INTEGER :: iyearprev  ! year of previous timestep 

    !INTEGER :: idayout         ! write out daily output
    !INTEGER :: imonthout! write out daily output
    !INTEGER :: iyearout 
    INTEGER :: nLndPts
    INTEGER :: i,k,ind,nint,IntSib,itr
    REAL (KIND=r8), PARAMETER   :: tice   =      271.16e0_r8! constant tice
    hc    =0.0_r8
    hg    =0.0_r8
    ec    =0.0_r8
    eg    =0.0_r8
    d     =0.0_r8
    bstar =0.0_r8

    !  idatec(1)....hour(00/12)
    !  idatec(2)....month
    !  idatec(3)....day of month
    !  idatec(4)....year
    iday      = idatec   (3) ! current day in month (1-31 or 1-30, or 1-28)
    imonth    = idatec   (2) ! current month (1 - 12)  
    iyear     = idatec   (4) ! current year 

    imonthprev= idateprev(2) ! month of previous timestep 
    !idayprev  = idateprev(3) ! day in month of previous timestep
    !iyearprev = idateprev(4) ! year of previous timestep 

    !calday = julday (imonth,iday,iyear,tod)
    rbyg  =rgas/grav*delsig*0.5_r8
    !CALL MsgOne(h," IBIS: IBIS Tendencies")
    !IF(myId ==0)WRITE(*,*)imonthprev
    nLndPts=0
    DO i=1,nCols
       GSW(I) = xvisb  (i)+xvisd  (i)+xnirb  (i)+xnird  (i)
       GLW(I) = fira2  (i)
       IF (iMask(i) >= 1_i8) THEN
          nLndPts=nLndPts+1 
          speedm (nLndPts  ) = MAX(SQRT((gu(i,1) /SIN(colrad(i)))**2  + (gv(i,1) /SIN(colrad(i)))**2),0.5_r8)
          ta     (nLndPts  ) = gt(i,1)
          qa     (nLndPts  ) = MAX(gq(i,1),0.00000012_r8) 
          ux     (nLndPts  ) = MAX(gu(i,1) /SIN(colrad(i)),0.5_r8)
          uy     (nLndPts  ) = MAX(gv(i,1) /SIN(colrad(i)),0.5_r8)
          psurf  (nLndPts  ) = gps    (i)*100.00_r8          ! surface pressure hPa -> Pa
          solad  (nLndPts,1) = xvisb  (i)  !! number of solar radiation wavebands : vis, nir
          solai  (nLndPts,1) = xvisd  (i)
          solad  (nLndPts,2) = xnirb  (i)
          solai  (nLndPts,2) = xnird  (i)
          fira   (nLndPts  ) = fira2  (i)      
          coszen (nLndPts  ) = zenith (i)
          raina  (nLndPts  ) = (ppli(i) + ppci(i) )*1.0e-3_r8  !convert mm/s to m/s
          snowa  (nLndPts  ) = snow   (i) *1.0e-3_r8  !convert mm/s to m/s
          ua     (nLndPts  ) = speedm (nLndPts) 
          q34    (nLndPts,jb) =  MAX(q34(nLndPts,jb),0.00000012_r8) 
          q12    (nLndPts,jb) =  MAX(q34(nLndPts,jb),0.00000012_r8) 
       END IF
    END DO
    InitMod = (initlz >= 0 .AND. ktm == -1 .AND. kt == 0 )
    IF(InitMod)THEN
       nLndPts=0
       DO i=1,nCols
          !HML  (i) = oml_hml0 + 2.0*SQRT((gu(i,1) /SIN(colrad(i)))**2  + (gv(i,1) /SIN(colrad(i)))**2) 
          IF( iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             tu0(nLndPts,jb) = ta     (nLndPts  )
             tu (nLndPts,jb) = ta     (nLndPts  )
             tum(nLndPts,jb) = ta     (nLndPts  )

             ts (nLndPts,jb) = ta     (nLndPts  )

             tl (nLndPts,jb) = ta     (nLndPts  )

             t12(nLndPts,jb) = ta     (nLndPts  )

             t34(nLndPts,jb) = ta     (nLndPts  )       

             tgm(nLndPts,jb) = tg     (nLndPts,jb) 

             q12(nLndPts,jb) = qa     (nLndPts  )       
             q34(nLndPts,jb) = qa     (nLndPts  )       
             DO k=1,nsoilay 
                tsoi0(nLndPts,k,jb) =ta     (nLndPts  )! tsoi (nLndPts,k,jb)
                tsoim(nLndPts,k,jb) =ta     (nLndPts  )! tsoi (nLndPts,k,jb)
                tsoi (nLndPts,k,jb) =ta     (nLndPts  )! tsoi (nLndPts,k,jb)
             END DO
             DO k=1,nsoilay
                wsoi0(nLndPts,k,jb)  = wsoi (nLndPts,k,jb)
                wsoim(nLndPts,k,jb)  = wsoi (nLndPts,k,jb)
                wsoi (nLndPts,k,jb)  = wsoi (nLndPts,k,jb)
             END DO
          END IF
       END DO
    END IF
    nLndPts=0
    DO i=1,nCols
       IF( iMask(i) >= 1_i8)THEN 
          nLndPts=nLndPts+1
          !          wliqmin  =0.00_r8         !local !
          !          wsnomin  =0.00_r8         !local !
          !          grunof   (nLndPts,jb)=0.00_r8 !local !
          !          ginvap   (nLndPts,jb)=0.00_r8 !local !
          !          gsuvap   (nLndPts,jb)=0.00_r8 !local !
          !          gtrans   (nLndPts,jb)=0.00_r8 !local !
          !          gtransu  (nLndPts,jb)=0.00_r8 !local !
          !          gtransl  (nLndPts,jb)=0.00_r8 !local !
          !          gdrain   (nLndPts,jb)=0.00_r8 !local !
          !          gadjust  (nLndPts,jb)=0.00_r8 !local !
          !          topparu  (nLndPts,jb)=0.00_r8 !local !
          !          topparl  (nLndPts,jb)=0.00_r8 !local !
          !          su       (nLndPts,jb)=0.00_r8 !local !
          !          ss       (nLndPts,jb)=0.00_r8 !local !
          !          sl       (nLndPts,jb)=0.00_r8 !local !
          !          agcub    (nLndPts,jb)=0.00_r8!local !
          !          agcuc    (nLndPts,jb)=0.00_r8!local !
          !          ancub    (nLndPts,jb)=0.00_r8!local !
          !          ancuc    (nLndPts,jb)=0.00_r8!local !
          !          totcondub(nLndPts,jb)=0.00_r8!local !
          !          totconduc(nLndPts,jb)=0.00_r8!local !
          !          agcls    (nLndPts,jb)=0.00_r8!local !
          !          agcl4    (nLndPts,jb)=0.00_r8!local !
          !          agcl3    (nLndPts,jb)=0.00_r8!local !
          !          ancls    (nLndPts,jb)=0.00_r8!local !
          !          ancl4    (nLndPts,jb)=0.00_r8!local !
          !          ancl3    (nLndPts,jb)=0.00_r8!local !
          !          totcondls(nLndPts,jb)=0.00_r8!local !
          !          totcondl3(nLndPts,jb)=0.00_r8!local !
          !         totcondl4(nLndPts,jb)=0.00_r8!local !
          !          consoi   (nLndPts,1:nsoilay,jb)=0.00_r8!local ! 
          !          qglif    (nLndPts,1:4,jb)=0.00_r8!local ! 
          !          hvasug   (nLndPts,jb)=0.00_r8!local ! 
          !          hvasui   (nLndPts,jb)=0.00_r8!local ! 
          !          stressl  (nLndPts,1:nsoilay,jb)=0.00_r8!local ! 
          !          stressu  (nLndPts,1:nsoilay,jb)=0.00_r8!local ! 
          !          stresstl (nLndPts,jb)=0.00_r8!local !  
          !          stresstu (nLndPts,jb)=0.00_r8!local !         
          !          upsoiu   (nLndPts,1:nsoilay,jb)=0.00_r8!local !         
          !          upsoil   (nLndPts,1:nsoilay,jb)=0.00_r8!local !         
          !          heatg    (nLndPts,jb)=0.00_r8!local !                  
          !          heati    (nLndPts,jb)=0.00_r8!local !                  
          !          asurd    (nLndPts,1:nband,jb)=0.00_r8!local !                 
          !          asuri    (nLndPts,1:nband,jb)=0.00_r8!local !                 
          !          tmm     (ncount)  = tmgm  (ncount,latco)
          !          tm0     (ncount)  = tmgp  (ncount,latco)
          DO k=1,nsoilay 
             tsoi0(nLndPts,k,jb) = tsoi0 (nLndPts,k,jb)
             tsoim(nLndPts,k,jb) = tsoim (nLndPts,k,jb)
             tsoi (nLndPts,k,jb) = tsoi  (nLndPts,k,jb)
          END DO
          DO k=1,nsoilay
             wsoi0(nLndPts,k,jb)  = wsoi0 (nLndPts,k,jb)
             wsoim(nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
             wsoi (nLndPts,k,jb)  = wsoi  (nLndPts,k,jb)
          END DO
       END IF
    END DO


    IF(InitMod)THEN
       nint=2
       IntSib=5
    ELSE
       nint=1
       IntSib=1
    END IF
    IF(TRIM(iswrad).NE.'NON'.AND.TRIM(ilwrad).NE.'NON') THEN
       IF(InitMod .AND. nLndPts >= 1)THEN
          DO ind=1,nint
             nLndPts=0
             DO i=1,nCols
                IF (iMask(i) >= 1_i8) THEN
                   nLndPts=nLndPts+1
                   !
                   !     precipitation
                   !
                   raina  (nLndPts  ) =0.0e0_r8*1.0e-3_r8  !convert mm/s to m/s
                END IF
             END DO
             DO itr=1,IntSib
                CALL IbisDrv (tod,pi,stef,vonk,grav,tmelt,hfus,hvap,hsub,ch2o,cice,cair,cvap,rair,rvap,cappa, &
                     rhow,nLndPts,nband,nsoilay,nsnolay,npft,epsilon,dtc3x,doalb,&
                     ginvap       (1:nLndPts,jb),gsuvap      (1:nLndPts,jb),gtrans         (1:nLndPts,jb), &
                     gtransu      (1:nLndPts,jb),gtransl        (1:nLndPts,jb),grunof         (1:nLndPts,jb),gdrain           (1:nLndPts,jb), &
                     gadjust      (1:nLndPts,jb),a10scalparamu(1:nLndPts,jb),a10daylightu(1:nLndPts,jb),a10scalparaml(1:nLndPts,jb), &
                     a10daylightl (1:nLndPts,jb),vmax_pft        (1:npft)    ,tau15                     ,kc15                       , &
                     ko15,cimax,gammaub,alpha3,theta3,beta3,coefmub,coefbub,gsubmin,gammauc,coefmuc,coefbuc, &
                     gsucmin,gammals,coefmls,coefbls,gslsmin,gammal3,coefml3,coefbl3, &
                     gsl3min,gammal4,alpha4,theta4,beta4,coefml4,coefbl4,gsl4min, &
                     wliqu        (1:nLndPts,jb),wliqumax                    ,wsnou         (1:nLndPts,jb),wsnoumax                 , &
                     tu              (1:nLndPts,jb),wliqs        (1:nLndPts,jb),wliqsmax                     ,wsnos           (1:nLndPts,jb), &
                     wsnosmax                  ,ts                (1:nLndPts,jb),wliql         (1:nLndPts,jb),wliqlmax                 , &  
                     wsnol        (1:nLndPts,jb),wsnolmax                    ,tl          (1:nLndPts,jb),topparu           (1:nLndPts,jb), &
                     topparl      (1:nLndPts,jb),fl                (1:nLndPts,jb),fu          (1:nLndPts,jb),lai      (1:nLndPts,1:2,jb), &
                     sai          (1:nLndPts,1:2,jb),rhoveg      (1:nband,1:2),tauveg        (1:nband,1:2),orieh    (1:2)               , &
                     oriev    (1:2)           ,wliqmin                    ,wsnomin                     ,t12           (1:nLndPts,jb), &
                     tdripu                   ,tblowu                    ,tdrips                     ,tblows                       , &
                     t34              (1:nLndPts,jb),tdripl                    ,tblowl                     ,ztop     (1:nLndPts,1:2,jb), &
                     alaiml                   ,zbot     (1:nLndPts,1:2,jb),alaimu                     ,froot        (1:nsoilay,1:2), &
                     q34              (1:nLndPts,jb),q12          (1:nLndPts,jb),su          (1:nLndPts,jb),cleaf                       , &
                     dleaf        (1:2)          ,ss                (1:nLndPts,jb),cstem                     ,dstem                       , &
                     sl              (1:nLndPts,jb),cgrass                    ,ciub         (1:nLndPts,jb),ciuc           (1:nLndPts,jb), &
                     exist (1:nLndPts,1:npft,jb),csub         (1:nLndPts,jb),gsub         (1:nLndPts,jb),csuc           (1:nLndPts,jb), &
                     gsuc              (1:nLndPts,jb),agcub        (1:nLndPts,jb),agcuc         (1:nLndPts,jb),ancub           (1:nLndPts,jb), &
                     ancuc        (1:nLndPts,jb),totcondub        (1:nLndPts,jb),totconduc   (1:nLndPts,jb),cils           (1:nLndPts,jb), &
                     cil3              (1:nLndPts,jb),cil4         (1:nLndPts,jb),csls         (1:nLndPts,jb),gsls           (1:nLndPts,jb), &
                     csl3              (1:nLndPts,jb),gsl3         (1:nLndPts,jb),csl4         (1:nLndPts,jb),gsl4           (1:nLndPts,jb), &
                     agcls        (1:nLndPts,jb),agcl4        (1:nLndPts,jb),agcl3         (1:nLndPts,jb),ancls           (1:nLndPts,jb), &
                     ancl4        (1:nLndPts,jb),ancl3        (1:nLndPts,jb),totcondls   (1:nLndPts,jb),totcondl3    (1:nLndPts,jb), &
                     totcondl4    (1:nLndPts,jb),chu,chs,chl,frac  (1:nLndPts,1:npft,jb),tlsub          (1:nLndPts,jb),z0sno,rhos, &
                     consno                   ,hsnotop                    ,hsnomin                     ,fimin                       , &
                     fimax                    ,fi                (1:nLndPts,jb),tsno(1:nLndPts,1:nsnolay,jb),hsno(1:nLndPts,1:nsnolay,jb),&
                     sand(1:nLndPts,1:nsoilay,jb),clay(1:nLndPts,1:nsoilay,jb),poros(1:nLndPts,1:nsoilay,jb),wsoi(1:nLndPts,1:nsoilay,jb), &
                     wisoi(1:nLndPts,1:nsoilay,jb),consoi(1:nLndPts,1:nsoilay,jb),zwpmax             ,wpud(1:nLndPts,jb)         , &
                     wipud        (1:nLndPts,jb),wpudmax                    ,qglif   (1:nLndPts,1:4,jb),tsoi(1:nLndPts,1:nsoilay,jb), &
                     hvasug       (1:nLndPts,jb),hvasui        (1:nLndPts,jb),albsav         (1:nLndPts,jb),albsan           (1:nLndPts,jb), &
                     tg              (1:nLndPts,jb),ti                (1:nLndPts,jb),z0soi(1:nLndPts,jb)       ,swilt(1:nLndPts,1:nsoilay,jb), &
                     sfield(1:nLndPts,1:nsoilay,jb),stressl(1:nLndPts,1:nsoilay,jb),stressu(1:nLndPts,1:nsoilay,jb),stresstl(1:nLndPts,jb), &
                     stresstu (1:nLndPts,jb)    ,csoi(1:nLndPts,1:nsoilay,jb),rhosoi(1:nLndPts,1:nsoilay,jb),hsoi(1:nsoilay+1)   , &
                     suction(1:nLndPts,1:nsoilay,jb),bex (1:nLndPts,1:nsoilay,jb),upsoiu(1:nLndPts,1:nsoilay,jb),upsoil(1:nLndPts,1:nsoilay,jb), &
                     heatg   (1:nLndPts,jb),heati   (1:nLndPts,jb),hydraul(1:nLndPts,1:nsoilay,jb),porosflo(1:nLndPts,1:nsoilay,jb), &
                     ibex(1:nLndPts,1:nsoilay,jb),bperm          ,hflo      (1:nLndPts,1:nsoilay+1,jb),ta       (1:nLndPts)       , &
                     asurd (1:nLndPts,1:nBand,jb),asuri(1:nLndPts,1:nBand,jb),coszen(1:nLndPts)              ,solad (1:nLndPts,1:nBand)  , &
                     solai (1:nLndPts,1:nBand)  ,fira      (1:nLndPts)      ,raina (1:nLndPts)             ,qa    (1:nLndPts)               , &
                     psurf    (1:nLndPts)          ,snowa     (1:nLndPts)      ,ua    (1:nLndPts)             ,o2conc                       , &
                     co2conc                  ,td         (1:nLndPts,jb)            ,vzero    (1:nLndPts,jb)   ,ndaypy                       , &
                     nppdummy(1:nLndPts,1:npft,jb) ,cbiow (1:nLndPts,1:npft,jb),sapfrac(1:nLndPts,jb)      ,cbior(1:nLndPts,1:npft,jb), &
                     tco2root(1:nLndPts,jb)          ,tneetot(1:nLndPts,jb)            ,tco2mic (1:nLndPts,jb)       ,a10td    (1:nLndPts,jb)    , &
                     a10ancub(1:nLndPts,jb)          ,a10ancuc(1:nLndPts,jb)     ,a10ancls(1:nLndPts,jb)    ,a10ancl3(1:nLndPts,jb)     , &
                     a10ancl4(1:nLndPts,jb)          ,ndtimes(1:nLndPts,jb), &
                     adrain  (1:nLndPts,jb)          ,adsnow  (1:nLndPts,jb),tnpptot(1:nLndPts,jb),adaet   (1:nLndPts,jb)    ,adtrunoff(1:nLndPts,jb)    , &
                     adsrunoff(1:nLndPts,jb)    ,addrainage(1:nLndPts,jb)   ,adrh    (1:nLndPts,jb)    ,adsnod   (1:nLndPts,jb)    , &
                     adsnof   (1:nLndPts,jb)    ,adwsoi    (1:nLndPts,jb)   ,adtsoi  (1:nLndPts,jb)    ,adwisoi  (1:nLndPts,jb)    , &
                     adtlaysoi(1:nLndPts,jb)    ,adwlaysoi (1:nLndPts,jb)   ,adwsoic (1:nLndPts,jb)    ,adtsoic  (1:nLndPts,jb)    , &
                     adco2mic (1:nLndPts,jb)    ,adco2root (1:nLndPts,jb)   ,adco2soi(1:nLndPts,jb)    ,adco2ratio(1:nLndPts,jb)   , &
                     adnmintot(1:nLndPts,jb)    ,decompl   (1:nLndPts,jb)   ,decomps (1:nLndPts,jb)    ,tnmin    (1:nLndPts,jb)    , &
                     ndaypm                     ,nmtimes(1:nLndPts,jb),amrain   (1:nLndPts,jb)    , &
                     amsnow    (1:nLndPts,jb)   ,amaet     (1:nLndPts,jb)   ,amtrunoff(1:nLndPts,jb)   ,amsrunoff(1:nLndPts,jb)    , &
                     amdrainage(1:nLndPts,jb)   ,amtemp    (1:nLndPts,jb)   ,amqa     (1:nLndPts,jb)    , &
                     amsolar   (1:nLndPts,jb)   ,amirup   (1:nLndPts,jb)   ,amirdown (1:nLndPts,jb)    , &
                     amsens    (1:nLndPts,jb)   ,amlatent  (1:nLndPts,jb)   ,amlaiu   (1:nLndPts,jb)   ,amlail   (1:nLndPts,jb)    , &
                     amtsoi    (1:nLndPts,jb)   ,amwsoi    (1:nLndPts,jb)   ,amwisoi  (1:nLndPts,jb)   ,amvwc    (1:nLndPts,jb)    , &
                     amawc     (1:nLndPts,jb)   ,amsnod    (1:nLndPts,jb)   ,amsnof   (1:nLndPts,jb)   ,amnpp(1:nLndPts,1:npft,jb) , &
                     amnpptot  (1:nLndPts,jb)   ,amco2mic  (1:nLndPts,jb)   ,amco2root(1:nLndPts,jb)   ,amco2soi (1:nLndPts,jb)    , &
                     amco2ratio(1:nLndPts,jb)   ,amneetot  (1:nLndPts,jb)   ,amnmintot(1:nLndPts,jb)   ,amalbedo (1:nLndPts,jb)    , &
                     amtsoil(1:nLndPts,1:nsoilay,jb) ,amwsoil(1:nLndPts,1:nsoilay,jb),amwisoil(1:nLndPts,1:nsoilay,jb), nytimes(1:nLndPts,jb), &
                     aysolar   (1:nLndPts,jb)   ,ayirup    (1:nLndPts,jb)   ,ayirdown (1:nLndPts,jb)   ,aysens        (1:nLndPts,jb)   , &
                     aylatent  (1:nLndPts,jb)   ,ayprcp    (1:nLndPts,jb)   ,ayaet    (1:nLndPts,jb)   ,aytrans        (1:nLndPts,jb)   , &
                     aytrunoff (1:nLndPts,jb)   ,aysrunoff (1:nLndPts,jb)   ,aydrainage(1:nLndPts,jb)  ,aydwtot        (1:nLndPts,jb)   , & 
                     aywsoi    (1:nLndPts,jb)   ,aywisoi   (1:nLndPts,jb)   ,aytsoi   (1:nLndPts,jb)   ,ayvwc        (1:nLndPts,jb)   , &
                     ayawc     (1:nLndPts,jb)   ,aystresstu(1:nLndPts,jb)   ,aystresstl(1:nLndPts,jb)  ,aygpp(1:nLndPts,1:npft,jb) , &
                     aygpptot  (1:nLndPts,jb)   ,aynpp(1:nLndPts,1:npft,jb) ,aynpptot (1:nLndPts,jb)   ,ayco2mic  (1:nLndPts,jb)   , &
                     ayco2root (1:nLndPts,jb)   ,ayco2soi  (1:nLndPts,jb)   ,ayneetot (1:nLndPts,jb)   ,ayrootbio (1:nLndPts,jb)   , &
                     aynmintot (1:nLndPts,jb)   ,ayalit    (1:nLndPts,jb)   ,ayblit   (1:nLndPts,jb)   ,aycsoi        (1:nLndPts,jb)   , & 
                     aycmic    (1:nLndPts,jb)   ,ayanlit   (1:nLndPts,jb)   ,aybnlit  (1:nLndPts,jb)   ,aynsoi        (1:nLndPts,jb)   , &
                     ayalbedo  (1:nLndPts,jb)   ,totalit   (1:nLndPts,jb)   ,totrlit  (1:nLndPts,jb)   ,totcsoi  (1:nLndPts,jb)    , &
                     totcmic   (1:nLndPts,jb)   , &
                     totanlit  (1:nLndPts,jb)   ,totrnlit  (1:nLndPts,jb)   ,totnsoi  (1:nLndPts,jb)   ,totnmic        (1:nLndPts,jb)   , &
                     totlit    (1:nLndPts,jb)   ,totfall   (1:nLndPts,jb)   ,totnlit  (1:nLndPts,jb)   ,firefac        (1:nLndPts,jb)   , &
                     wtot      (1:nLndPts,jb)   ,storedn   (1:nLndPts,jb)   ,yrleach  (1:nLndPts,jb)   ,ynleach        (1:nLndPts,jb)   , & 
                     falll     (1:nLndPts,jb)   ,fallr     (1:nLndPts,jb)   ,fallw    (1:nLndPts,jb)   ,clitlm        (1:nLndPts,jb)   , &
                     clitls    (1:nLndPts,jb)   ,clitrm    (1:nLndPts,jb)   ,clitrs   (1:nLndPts,jb)   ,clitwm        (1:nLndPts,jb)   , &
                     clitws    (1:nLndPts,jb)   ,csoislop  (1:nLndPts,jb)   ,csoislon (1:nLndPts,jb)   ,csoipas        (1:nLndPts,jb)   , &
                     clitll    (1:nLndPts,jb)   ,clitrl    (1:nLndPts,jb)   ,clitwl   (1:nLndPts,jb)   ,tc        (1:nLndPts,jb)   , &
                     agddu     (1:nLndPts,jb)   ,tempu     (1:nLndPts,jb)   ,agddl    (1:nLndPts,jb)   ,templ        (1:nLndPts,jb)   , &
                     dropu     (1:nLndPts,jb)   ,dropls    (1:nLndPts,jb)   ,dropl4   (1:nLndPts,jb)   ,dropl3        (1:nLndPts,jb)   , &
                     plai(1:nLndPts,1:npft,jb)  ,iday,imonth,iyear,iyear0,isimveg,spinmax, &
                     amts2     (1:nLndPts,jb)   ,amtransu  (1:nLndPts,jb)  ,amtransl (1:nLndPts,jb),   amsuvap  (1:nLndPts,jb)        ,&
                     aminvap   (1:nLndPts,jb)   ,ux             (1:nLndPts)      ,uy       (1:nLndPts)      ,taux        (1:nLndPts)      , &
                     tauy           (1:nLndPts),ts2(1:nLndPts),qs2(1:nLndPts),deltat   (1:nLndPts,jb)   ,gdd0        (1:nLndPts,jb)   , &  
                     gdd0this  (1:nLndPts,jb)   ,tcthis    (1:nLndPts,jb)   ,twthis   (1:nLndPts,jb)   ,tcmin        (1:nLndPts,jb)   , &
                     gdd5           (1:nLndPts,jb)   ,gdd5this  (1:nLndPts,jb)   ,TminL    (1:npft)       ,TminU        (1:npft)       , &
                     Twarm     (1:npft)          ,GDD       (1:npft)            ,aleaf    (1:npft)       ,awood        (1:npft)       , &
                     cbiol(1:nLndPts,1:npft,jb) ,aroot     (1:npft)            ,disturbf (1:nLndPts)      ,disturbo  (1:nLndPts)      , &
                     specla    (1:npft)          ,biomass(1:nLndPts,1:npft,jb),totlaiu (1:nLndPts,jb)   ,totlail        (1:nLndPts,jb)   , &  
                     totbiou   (1:nLndPts,jb)   ,totbiol   (1:nLndPts,jb)   ,woodnorm                     ,vegtype0  (1:nLndPts,jb)   , &
                     tauwood0  (1:npft)          ,tauwood   (1:nLndPts,1:npft,jb)            ,tauleaf  (1:npft)       ,tauroot        (1:npft)       , &
                     xminlai,cdisturb  (1:nLndPts,jb)   ,ayanpp(1:nLndPts,1:npft,jb),ayanpptot(1:nLndPts,jb),jdt,&
                     imonthprev,isimco2,isimfire, co2init,tw(1:nLndPts,jb),&
                     fvapa  (1:nLndPts),fsena (1:nLndPts),z0(1:nLndPts),ustar(1:nLndPts), hc(1:nLndPts),hg (1:nLndPts),ec (1:nLndPts),&
                     eg (1:nLndPts) ,d (1:nLndPts) ,cu (1:nLndPts),firb  (1:nLndPts),tgpptot (1:nLndPts),bstar1(1:nLndPts))

                nLndPts=0
                DO i=1,nCols
                   IF (iMask(i) >= 1_i8) THEN
                      nLndPts=nLndPts+1
                      tu (nLndPts,jb) = ta     (nLndPts  )
                      ts (nLndPts,jb) = ta     (nLndPts  )
                      tl (nLndPts,jb) = ta     (nLndPts  )
                      t12(nLndPts,jb) = ta     (nLndPts  )
                      t34(nLndPts,jb) = ta     (nLndPts  )        
                      tgm(nLndPts,jb) = tg     (nLndPts,jb) 
                      q12(nLndPts,jb) = qa     (nLndPts  )        
                      q34(nLndPts,jb) = qa     (nLndPts  )        
                   END IF
                END DO

             END DO
             DO k=1,nsoilay
                DO i=1,nLndPts
                   tsoim(i,k,jb) = tsoi (i,k,jb)
                   tsoi (i,k,jb) = tsoim(i,k,jb)
                END DO
             END DO
             DO k=1,nsoilay
                DO i=1,nLndPts
                   wsoim(i,k,jb) = wsoi(i,k,jb)
                   wsoi(i,k,jb)  = wsoim(i,k,jb)
                END DO
             END DO
             DO i=1,nLndPts
                !    capac(i,1)=capacm(i,1) wliqu
                !    capac(i,2)=capacm(i,2) wliqu
                !    tu   (i,jb)  =tum   (i,jb)
                IF(ind.EQ.1) THEN
                   tmin (i) =tg (i,jb)
                ELSE
                   tmax (i) =tg (i,jb)
                END IF
                tgm   (i,jb)=tg (i,jb)
                tg   (i,jb) =tgm(i,jb)
             END DO
          END DO
          DO k=1,nsoilay
             DO i=1,nLndPts
                !          td   (i,k) =tdm   (i,k)
                tsoi (i,k,jb) = tsoim(i,k,jb)
                !    td   (i,k) =0.9_r8*0.5_r8*(tmax(i)+tmin(i))+0.1_r8*tdm(i,k)
                tsoi (i,k,jb) = 0.9_r8*0.5_r8*(tmax(i)+tmin(i))+0.1_r8*tsoim(i,k,jb)
                !    tdm  (i,k) =td(i,k)
                tsoim(i,k,jb) = tsoi (i,k,jb)
                !    td0  (i,k) =td(i,k)
                tsoi0(i,k,jb) = tsoi (i,k,jb)
             END DO
          END DO
          !
          !     this is a start of equilibrium tg,tc comp.
          !
          DO i=1,nLndPts
             IF(coszen(i).LT.0.0e0_r8) THEN
                tgm  (i,jb)  =tmin(i)
                tg0  (i,jb)  =tmin(i)
             END IF
          END DO
       END IF
    END IF
    !PPPPP
    IF(nLndPts.GE.1) THEN
       nLndPts=0
       DO i=1,nCols
          IF (iMask(i) >= 1_i8) THEN
             nLndPts=nLndPts+1
             !
             !     precipitation
             !
             raina  (nLndPts  ) = (ppli(i) + ppci(i) )*1.0e-3_r8  !convert mm/s to m/s
          END IF
       END DO
       CALL IbisDrv (tod,pi,stef,vonk,grav,tmelt,hfus,hvap,hsub,ch2o,cice,cair,cvap,rair,rvap,cappa, &
            rhow,nLndPts,nband,nsoilay,nsnolay,npft,epsilon,dtc3x,doalb,&
            ginvap       (1:nLndPts,jb),gsuvap      (1:nLndPts,jb),gtrans         (1:nLndPts,jb), &
            gtransu      (1:nLndPts,jb),gtransl        (1:nLndPts,jb),grunof         (1:nLndPts,jb),gdrain           (1:nLndPts,jb), &
            gadjust      (1:nLndPts,jb),a10scalparamu(1:nLndPts,jb),a10daylightu(1:nLndPts,jb),a10scalparaml(1:nLndPts,jb), &
            a10daylightl (1:nLndPts,jb),vmax_pft        (1:npft)    ,tau15                     ,kc15                       , &
            ko15,cimax,gammaub,alpha3,theta3,beta3,coefmub,coefbub,gsubmin,gammauc,coefmuc,coefbuc, &
            gsucmin,gammals,coefmls,coefbls,gslsmin,gammal3,coefml3,coefbl3, &
            gsl3min,gammal4,alpha4,theta4,beta4,coefml4,coefbl4,gsl4min, &
            wliqu        (1:nLndPts,jb),wliqumax                    ,wsnou         (1:nLndPts,jb),wsnoumax                 , &
            tu              (1:nLndPts,jb),wliqs        (1:nLndPts,jb),wliqsmax                     ,wsnos           (1:nLndPts,jb), &
            wsnosmax                  ,ts                (1:nLndPts,jb),wliql         (1:nLndPts,jb),wliqlmax                 , &  
            wsnol        (1:nLndPts,jb),wsnolmax                    ,tl          (1:nLndPts,jb),topparu           (1:nLndPts,jb), &
            topparl      (1:nLndPts,jb),fl                (1:nLndPts,jb),fu          (1:nLndPts,jb),lai      (1:nLndPts,1:2,jb), &
            sai          (1:nLndPts,1:2,jb),rhoveg      (1:nband,1:2),tauveg        (1:nband,1:2),orieh    (1:2)               , &
            oriev    (1:2)           ,wliqmin                    ,wsnomin                     ,t12           (1:nLndPts,jb), &
            tdripu                   ,tblowu                    ,tdrips                     ,tblows                       , &
            t34              (1:nLndPts,jb),tdripl                    ,tblowl                     ,ztop     (1:nLndPts,1:2,jb), &
            alaiml                   ,zbot     (1:nLndPts,1:2,jb),alaimu                     ,froot        (1:nsoilay,1:2), &
            q34              (1:nLndPts,jb),q12          (1:nLndPts,jb),su          (1:nLndPts,jb),cleaf                       , &
            dleaf        (1:2)          ,ss                (1:nLndPts,jb),cstem                     ,dstem                       , &
            sl              (1:nLndPts,jb),cgrass                    ,ciub         (1:nLndPts,jb),ciuc           (1:nLndPts,jb), &
            exist (1:nLndPts,1:npft,jb),csub         (1:nLndPts,jb),gsub         (1:nLndPts,jb),csuc           (1:nLndPts,jb), &
            gsuc              (1:nLndPts,jb),agcub        (1:nLndPts,jb),agcuc         (1:nLndPts,jb),ancub           (1:nLndPts,jb), &
            ancuc        (1:nLndPts,jb),totcondub        (1:nLndPts,jb),totconduc   (1:nLndPts,jb),cils           (1:nLndPts,jb), &
            cil3              (1:nLndPts,jb),cil4         (1:nLndPts,jb),csls         (1:nLndPts,jb),gsls           (1:nLndPts,jb), &
            csl3              (1:nLndPts,jb),gsl3         (1:nLndPts,jb),csl4         (1:nLndPts,jb),gsl4           (1:nLndPts,jb), &
            agcls        (1:nLndPts,jb),agcl4        (1:nLndPts,jb),agcl3         (1:nLndPts,jb),ancls           (1:nLndPts,jb), &
            ancl4        (1:nLndPts,jb),ancl3        (1:nLndPts,jb),totcondls   (1:nLndPts,jb),totcondl3    (1:nLndPts,jb), &
            totcondl4    (1:nLndPts,jb),chu,chs,chl,frac  (1:nLndPts,1:npft,jb),tlsub          (1:nLndPts,jb),z0sno,rhos, &
            consno                   ,hsnotop                    ,hsnomin                     ,fimin                       , &
            fimax                    ,fi                (1:nLndPts,jb),tsno(1:nLndPts,1:nsnolay,jb),hsno(1:nLndPts,1:nsnolay,jb),&
            sand(1:nLndPts,1:nsoilay,jb),clay(1:nLndPts,1:nsoilay,jb),poros(1:nLndPts,1:nsoilay,jb),wsoi(1:nLndPts,1:nsoilay,jb), &
            wisoi(1:nLndPts,1:nsoilay,jb),consoi(1:nLndPts,1:nsoilay,jb),zwpmax             ,wpud(1:nLndPts,jb)         , &
            wipud        (1:nLndPts,jb),wpudmax                    ,qglif   (1:nLndPts,1:4,jb),tsoi(1:nLndPts,1:nsoilay,jb), &
            hvasug       (1:nLndPts,jb),hvasui        (1:nLndPts,jb),albsav         (1:nLndPts,jb),albsan           (1:nLndPts,jb), &
            tg              (1:nLndPts,jb),ti                (1:nLndPts,jb),z0soi(1:nLndPts,jb)       ,swilt(1:nLndPts,1:nsoilay,jb), &
            sfield(1:nLndPts,1:nsoilay,jb),stressl(1:nLndPts,1:nsoilay,jb),stressu(1:nLndPts,1:nsoilay,jb),stresstl(1:nLndPts,jb), &
            stresstu (1:nLndPts,jb)    ,csoi(1:nLndPts,1:nsoilay,jb),rhosoi(1:nLndPts,1:nsoilay,jb),hsoi(1:nsoilay+1)   , &
            suction(1:nLndPts,1:nsoilay,jb),bex (1:nLndPts,1:nsoilay,jb),upsoiu(1:nLndPts,1:nsoilay,jb),upsoil(1:nLndPts,1:nsoilay,jb), &
            heatg   (1:nLndPts,jb),heati   (1:nLndPts,jb),hydraul(1:nLndPts,1:nsoilay,jb),porosflo(1:nLndPts,1:nsoilay,jb), &
            ibex(1:nLndPts,1:nsoilay,jb),bperm          ,hflo      (1:nLndPts,1:nsoilay+1,jb),ta       (1:nLndPts)       , &
            asurd (1:nLndPts,1:nBand,jb),asuri(1:nLndPts,1:nBand,jb),coszen(1:nLndPts)              ,solad (1:nLndPts,1:nBand)  , &
            solai (1:nLndPts,1:nBand)  ,fira      (1:nLndPts)      ,raina (1:nLndPts)             ,qa    (1:nLndPts)               , &
            psurf    (1:nLndPts)          ,snowa     (1:nLndPts)      ,ua    (1:nLndPts)             ,o2conc                       , &
            co2conc                  ,td         (1:nLndPts,jb)            ,vzero    (1:nLndPts,jb)   ,ndaypy                       , &
            nppdummy(1:nLndPts,1:npft,jb) ,cbiow (1:nLndPts,1:npft,jb),sapfrac(1:nLndPts,jb)      ,cbior(1:nLndPts,1:npft,jb), &
            tco2root(1:nLndPts,jb)          ,tneetot(1:nLndPts,jb)            ,tco2mic (1:nLndPts,jb)       ,a10td    (1:nLndPts,jb)    , &
            a10ancub(1:nLndPts,jb)          ,a10ancuc(1:nLndPts,jb)     ,a10ancls(1:nLndPts,jb)    ,a10ancl3(1:nLndPts,jb)     , &
            a10ancl4(1:nLndPts,jb)          ,ndtimes(1:nLndPts,jb), &
            adrain  (1:nLndPts,jb)          ,adsnow  (1:nLndPts,jb),tnpptot(1:nLndPts,jb),adaet   (1:nLndPts,jb)    ,adtrunoff(1:nLndPts,jb), &
            adsrunoff(1:nLndPts,jb)    ,addrainage(1:nLndPts,jb)   ,adrh    (1:nLndPts,jb)    ,adsnod   (1:nLndPts,jb)    , &
            adsnof   (1:nLndPts,jb)    ,adwsoi    (1:nLndPts,jb)   ,adtsoi  (1:nLndPts,jb)    ,adwisoi  (1:nLndPts,jb)    , &
            adtlaysoi(1:nLndPts,jb)    ,adwlaysoi (1:nLndPts,jb)   ,adwsoic (1:nLndPts,jb)    ,adtsoic  (1:nLndPts,jb)    , &
            adco2mic (1:nLndPts,jb)    ,adco2root (1:nLndPts,jb)   ,adco2soi(1:nLndPts,jb)    ,adco2ratio(1:nLndPts,jb)   , &
            adnmintot(1:nLndPts,jb)    ,decompl   (1:nLndPts,jb)   ,decomps (1:nLndPts,jb)    ,tnmin    (1:nLndPts,jb)    , &
            ndaypm                   ,nmtimes(1:nLndPts,jb),amrain   (1:nLndPts,jb)    , &
            amsnow    (1:nLndPts,jb)   ,amaet     (1:nLndPts,jb)   ,amtrunoff(1:nLndPts,jb)   ,amsrunoff(1:nLndPts,jb)    , &
            amdrainage(1:nLndPts,jb)   ,amtemp    (1:nLndPts,jb)   ,amqa     (1:nLndPts,jb)    , &
            amsolar   (1:nLndPts,jb)   ,amirup   (1:nLndPts,jb)   ,amirdown (1:nLndPts,jb)    , &
            amsens    (1:nLndPts,jb)   ,amlatent  (1:nLndPts,jb)   ,amlaiu   (1:nLndPts,jb)   ,amlail   (1:nLndPts,jb)    , &
            amtsoi    (1:nLndPts,jb)   ,amwsoi    (1:nLndPts,jb)   ,amwisoi  (1:nLndPts,jb)   ,amvwc    (1:nLndPts,jb)    , &
            amawc     (1:nLndPts,jb)   ,amsnod    (1:nLndPts,jb)   ,amsnof   (1:nLndPts,jb)   ,amnpp(1:nLndPts,1:npft,jb) , &
            amnpptot  (1:nLndPts,jb)   ,amco2mic  (1:nLndPts,jb)   ,amco2root(1:nLndPts,jb)   ,amco2soi (1:nLndPts,jb)    , &
            amco2ratio(1:nLndPts,jb)   ,amneetot  (1:nLndPts,jb)   ,amnmintot(1:nLndPts,jb)   ,amalbedo (1:nLndPts,jb)    , &
            amtsoil(1:nLndPts,1:nsoilay,jb) ,amwsoil(1:nLndPts,1:nsoilay,jb),amwisoil(1:nLndPts,1:nsoilay,jb), nytimes(1:nLndPts,jb), &
            aysolar   (1:nLndPts,jb)   ,ayirup    (1:nLndPts,jb)   ,ayirdown (1:nLndPts,jb)   ,aysens        (1:nLndPts,jb)   , &
            aylatent  (1:nLndPts,jb)   ,ayprcp    (1:nLndPts,jb)   ,ayaet    (1:nLndPts,jb)   ,aytrans        (1:nLndPts,jb)   , &
            aytrunoff (1:nLndPts,jb)   ,aysrunoff (1:nLndPts,jb)   ,aydrainage(1:nLndPts,jb)  ,aydwtot        (1:nLndPts,jb)   , & 
            aywsoi    (1:nLndPts,jb)   ,aywisoi   (1:nLndPts,jb)   ,aytsoi   (1:nLndPts,jb)   ,ayvwc        (1:nLndPts,jb)   , &
            ayawc     (1:nLndPts,jb)   ,aystresstu(1:nLndPts,jb)   ,aystresstl(1:nLndPts,jb)  ,aygpp(1:nLndPts,1:npft,jb) , &
            aygpptot  (1:nLndPts,jb)   ,aynpp(1:nLndPts,1:npft,jb) ,aynpptot (1:nLndPts,jb)   ,ayco2mic  (1:nLndPts,jb)   , &
            ayco2root (1:nLndPts,jb)   ,ayco2soi  (1:nLndPts,jb)   ,ayneetot (1:nLndPts,jb)   ,ayrootbio (1:nLndPts,jb)   , &
            aynmintot (1:nLndPts,jb)   ,ayalit    (1:nLndPts,jb)   ,ayblit   (1:nLndPts,jb)   ,aycsoi        (1:nLndPts,jb)   , & 
            aycmic    (1:nLndPts,jb)   ,ayanlit   (1:nLndPts,jb)   ,aybnlit  (1:nLndPts,jb)   ,aynsoi        (1:nLndPts,jb)   , &
            ayalbedo  (1:nLndPts,jb)   ,totalit   (1:nLndPts,jb)   ,totrlit  (1:nLndPts,jb)   ,totcsoi  (1:nLndPts,jb)    , &
            totcmic   (1:nLndPts,jb)   , &
            totanlit  (1:nLndPts,jb)   ,totrnlit  (1:nLndPts,jb)   ,totnsoi  (1:nLndPts,jb)   ,totnmic        (1:nLndPts,jb)   , &
            totlit    (1:nLndPts,jb)   ,totfall   (1:nLndPts,jb)   ,totnlit  (1:nLndPts,jb)   ,firefac        (1:nLndPts,jb)   , &
            wtot           (1:nLndPts,jb)   ,storedn   (1:nLndPts,jb)   ,yrleach  (1:nLndPts,jb)   ,ynleach        (1:nLndPts,jb)   , & 
            falll     (1:nLndPts,jb)   ,fallr     (1:nLndPts,jb)   ,fallw    (1:nLndPts,jb)   ,clitlm        (1:nLndPts,jb)   , &
            clitls    (1:nLndPts,jb)   ,clitrm    (1:nLndPts,jb)   ,clitrs   (1:nLndPts,jb)   ,clitwm        (1:nLndPts,jb)   , &
            clitws    (1:nLndPts,jb)   ,csoislop  (1:nLndPts,jb)   ,csoislon (1:nLndPts,jb)   ,csoipas        (1:nLndPts,jb)   , &
            clitll    (1:nLndPts,jb)   ,clitrl    (1:nLndPts,jb)   ,clitwl   (1:nLndPts,jb)   ,tc        (1:nLndPts,jb)   , &
            agddu     (1:nLndPts,jb)   ,tempu     (1:nLndPts,jb)   ,agddl    (1:nLndPts,jb)   ,templ        (1:nLndPts,jb)   , &
            dropu     (1:nLndPts,jb)   ,dropls    (1:nLndPts,jb)   ,dropl4   (1:nLndPts,jb)   ,dropl3        (1:nLndPts,jb)   , &
            plai(1:nLndPts,1:npft,jb)  ,iday,imonth,iyear,iyear0,isimveg,spinmax, &
            amts2     (1:nLndPts,jb)   ,amtransu  (1:nLndPts,jb)  ,amtransl (1:nLndPts,jb),   amsuvap  (1:nLndPts,jb)        ,&
            aminvap   (1:nLndPts,jb)   ,ux             (1:nLndPts)      ,uy       (1:nLndPts)      ,taux        (1:nLndPts)      , &
            tauy(1:nLndPts),ts2(1:nLndPts),qs2(1:nLndPts),deltat(1:nLndPts,jb),gdd0(1:nLndPts,jb), &  
            gdd0this  (1:nLndPts,jb)   ,tcthis    (1:nLndPts,jb)   ,twthis   (1:nLndPts,jb)   ,tcmin        (1:nLndPts,jb)   , &
            gdd5           (1:nLndPts,jb)   ,gdd5this  (1:nLndPts,jb)   ,TminL    (1:npft)       ,TminU        (1:npft)       , &
            Twarm     (1:npft)          ,GDD       (1:npft)            ,aleaf    (1:npft)       ,awood        (1:npft)       , &
            cbiol(1:nLndPts,1:npft,jb) ,aroot     (1:npft)            ,disturbf (1:nLndPts)      ,disturbo  (1:nLndPts)      , &
            specla    (1:npft)          ,biomass(1:nLndPts,1:npft,jb),totlaiu (1:nLndPts,jb)   ,totlail        (1:nLndPts,jb)   , &  
            totbiou   (1:nLndPts,jb)   ,totbiol   (1:nLndPts,jb)   ,woodnorm                     ,vegtype0  (1:nLndPts,jb)   , &
            tauwood0  (1:npft)          ,tauwood   (1:nLndPts,1:npft,jb)            ,tauleaf  (1:npft)       ,tauroot        (1:npft)       , &
            xminlai                  ,cdisturb  (1:nLndPts,jb)   ,ayanpp(1:nLndPts,1:npft,jb),ayanpptot(1:nLndPts,jb)   , &
            jdt,imonthprev,isimco2,isimfire,&
            co2init,tw(1:nLndPts,jb),fvapa  (1:nLndPts),fsena (1:nLndPts),z0(1:nLndPts),ustar(1:nLndPts) ,hc(1:nLndPts),&
            hg (1:nLndPts),ec (1:nLndPts),eg (1:nLndPts),d(1:nLndPts),cu(1:nLndPts),firb  (1:nLndPts),tgpptot (1:nLndPts),bstar1(1:nLndPts))


    END IF

    !
    !     sib time integaration and time filter
    !
    DO i=1,nLndPts
       !      qm(i)=MAX(1.0e-12_r8,qm(i))
    END DO

    nLndPts=0
    DO i=1,nCols
       IF (iMask(i) >= 1_i8) THEN
          nLndPts=nLndPts+1
          !tgeff(i)=SQRT ( SQRT (( firb(i)*stef )))fu fi fl        tgeff(i)=SQRT ( SQRT (( firb(i)*stef )))
          tsea(i)            =SQRT ( SQRT (( firb(nLndPts)*(1.0_r8  /stef) )))!(tu(nLndPts,jb) + tl(nLndPts,jb)+ti(nLndPts,jb) + tg(nLndPts,jb))/4.0_r8
          IF(omlmodel)TSK (i)=SQRT ( SQRT (( firb(nLndPts)*(1.0_r8  /stef) )))!(tu(nLndPts,jb) + tl(nLndPts,jb)+ti(nLndPts,jb) + tg(nLndPts,jb))/4.0_r8
       END IF
    END DO
    nLndPts=0
    DO i=1,nCols
       IF (iMask(i) >= 1_i8)THEN 
          nLndPts=nLndPts+1
          sens    (i) =  -fsena (nLndPts)
          evap    (i) =  -fvapa (nLndPts)*hvap

          tgm  (nLndPts,jb  ) = tgm     (nLndPts,jb )
          tg   (nLndPts,jb  ) = tg      (nLndPts,jb )
          tg0  (nLndPts,jb  ) = tg0     (nLndPts,jb )

          DO k=1,nsoilay 
             tsoi0(nLndPts,k,jb) = tsoi0 (nLndPts,k,jb)
             tsoim(nLndPts,k,jb) = tsoim (nLndPts,k,jb)
             tsoi (nLndPts,k,jb) = tsoi  (nLndPts,k,jb)
          END DO
          DO k=1,nsoilay
             wsoi0(nLndPts,k,jb)  = wsoi0 (nLndPts,k,jb)
             wsoim(nLndPts,k,jb)  = wsoim (nLndPts,k,jb)
             wsoi (nLndPts,k,jb)  = wsoi  (nLndPts,k,jb)
          END DO
       END IF
    END DO
    !
    !     sea or sea ice
    ! gu gv gps colrad sigki delsig sens evap umom vmom rmi rhi cond stor zorl rnet ztn2 THETA_2M VELC_2m MIXQ_2M
    ! THETA_10M VELC_10M MIXQ_10M
    ! mmax=ncols-nmax+1
    ! including case 1D physics

    IF(initlz.GE.0.AND.kt.EQ.0.AND.jdt.EQ.1)THEN
       Tsfc0(1:nCols,jb)=gt  (1:nCols,1)
       Qsfc0(1:nCols,jb)=gq  (1:nCols,1)
       Tsfcm(1:nCols,jb)=gt  (1:nCols,1)
       Qsfcm(1:nCols,jb)=gq  (1:nCols,1)
       tsfc (1:nCols)=gt  (1:nCols,1)
       qsfc (1:nCols)=gq  (1:nCols,1)
    END IF
    DO i=1,nCols
       IF(mskant(i) == 1_i8)THEN
          xsea (i) = tseam(i)
          tsfc (i) = Tsfcm(i,jb)
          qsfc (i) = Qsfcm(i,jb)
       END IF
    END DO

    IF(TRIM(OCFLUX) == 'WGFS')THEN
    CALL seasfc_wgfs( tmtx (1:nCols,1:3)  ,umtx (1:nCols,1:4)  ,qmtx(1:nCols,1:3)    ,& 
         slrad       (1:nCols)            ,&
         tsurf (1:nCols)        ,qsurf (1:nCols)           ,gu    (1:nCols,1:kMax)   ,&
         gv    (1:nCols,1:kMax) ,gt    (1:nCols,1:kMax)    ,gq    (1:nCols,1:kMax)   ,&
         gps   (1:nCols)        ,xsea  (1:nCols)           ,dtc3x                    ,&
         SIN(colrad(1:nCols))   ,sigki                     ,delsig                   ,&
         sens  (1:nCols)        ,evap  (1:nCols)           ,umom  (1:nCols)          ,&
         vmom  (1:nCols)        ,rmi   (1:nCols)           ,rhi   (1:nCols)          ,&
         cond  (1:nCols)        ,stor  (1:nCols)           ,zorl  (1:nCols)          ,&
         rnet  (1:nCols)        ,nCols                     ,kMax                     ,&
         Ustarm(1:nCols)        ,z0sea (1:nCols)           ,rho   (1:nCols)          ,&
         qsfc  (1:nCols)        ,tsfc  (1:nCols)           ,mskant(1:nCols)          ,&
         bstar (1:nCols)        ,iMask (1:nCols)           ,HML   (1:nCols)          ,&
         HUML  (1:nCols)        ,HVML  (1:nCols)           ,TSK   (1:nCols)          ,&
         GSW   (1:nCols)        ,GLW   (1:nCols)             )
    ELSE IF (TRIM(OCFLUX) == 'UKME')THEN
    CALL seasfc_ukme(tmtx (1:nCols,1:3)  ,umtx (1:nCols,1:4)  ,qmtx(1:nCols,1:3)   ,& 
         slrad       (1:nCols)            ,&
         tsurf (1:nCols)        ,qsurf (1:nCols)    ,gu    (1:nCols,1:kMax)          ,&
         gv    (1:nCols,1:kMax)        ,gt    (1:nCols,1:kMax)    ,gq    (1:nCols,1:kMax)          ,&
         gps   (1:nCols)        ,xsea  (1:nCols)    ,dtc3x                    ,&
         SIN(colrad(1:nCols))   ,sigki              ,delsig                   ,&
         sens  (1:nCols)        ,evap  (1:nCols)    ,umom  (1:nCols)          ,&
         vmom  (1:nCols)        ,rmi   (1:nCols)    ,rhi   (1:nCols)          ,&
         cond  (1:nCols)        ,stor  (1:nCols)    ,zorl  (1:nCols)          ,&
         rnet  (1:nCols)        ,nCols              ,kMax                     ,&
         Ustarm(1:nCols)        ,z0sea (1:nCols)    ,rho   (1:nCols)          ,&
         qsfc  (1:nCols)        ,tsfc  (1:nCols)    ,mskant(1:nCols)          ,&
         bstar (1:nCols)        ,iMask (1:nCols)    ,cldtot  (1:nCols,1:kMax) ,&
         ySwSfcNet(1:nCols)     ,LwSfcNet(1:nCols)  ,pblh(1:nCols)            ,&
         QCF (1:nCols,1:kMax)   ,QCL  (1:nCols,1:kMax) )
    ELSE IF (TRIM(OCFLUX) == 'COLA')THEN
    CALL seasfc_cola(tmtx (1:nCols,1:3)  ,umtx (1:nCols,1:4)  ,qmtx(1:nCols,1:3)   ,& 
         slrad       (1:nCols)            ,&
         tsurf (1:nCols)        ,qsurf (1:nCols)    ,gu    (1:nCols,1:kMax)          ,&
         gv    (1:nCols,1:kMax)        ,gt    (1:nCols,1:kMax)    ,gq    (1:nCols,1:kMax)          ,&
         gps   (1:nCols)        ,xsea  (1:nCols)    ,dtc3x                    ,&
         SIN(colrad(1:nCols))   ,sigki              ,delsig                   ,&
         sens  (1:nCols)        ,evap  (1:nCols)    ,umom  (1:nCols)          ,&
         vmom  (1:nCols)        ,rmi   (1:nCols)    ,rhi   (1:nCols)          ,&
         cond  (1:nCols)        ,stor  (1:nCols)    ,zorl  (1:nCols)          ,&
         rnet  (1:nCols)        ,nCols              ,kMax                     ,&
         Ustarm(1:nCols)        ,z0sea (1:nCols)    ,rho   (1:nCols)          ,&
         qsfc  (1:nCols)        ,tsfc  (1:nCols)    ,mskant(1:nCols)          ,&
         bstar (1:nCols)         )
    ELSE
       WRITE(*,*)'ERRO seasfc',OCFLUX
    END IF
    
    DO i=1,nCols
       IF(mskant(i) == 1_i8 .AND. tsea(i) <= 0.0e0_r8 .AND. tsurf(i) < tice+0.01e0_r8 ) THEN
          IF(intg.EQ.2) THEN
             IF(istrt.EQ.0) THEN
                tseam(i)=filta*tsea (i) + epsflt*(tseam(i)+xsea(i))
                qsfc (i)=MAX(1.0e-12_r8,qsfc(i))
                Tsfcm(i,jb)=filta*Tsfc0 (i,jb) + epsflt*(Tsfcm(i,jb)+tsfc(i))
                Qsfcm(i,jb)=filta*Qsfc0 (i,jb) + epsflt*(Qsfcm(i,jb)+qsfc(i))
             END IF
             tsea (i) = xsea(i)
             qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
             Tsfc0(i,jb) = tsfc(i)
             Qsfc0(i,jb) = qsfc(i)
          ELSE
             tsea (i) = xsea(i)
             tseam(i) = xsea(i)
             qsfc (i) = MAX(1.0e-12_r8,qsfc(i))
             Tsfc0(i,jb) = tsfc(i)
             Qsfc0(i,jb) = qsfc(i)
             Tsfcm(i,jb) = tsfc(i)
             Qsfcm(i,jb) = qsfc(i)
          END IF
       END IF
       IF(mskant(i) == 1_i8 .AND. tsea(i).LT.0.0e0_r8.AND.tsurf(i).GE.tice+0.01e0_r8) THEN
          tseam(i) = tsea (i)
          Tsfcm(i,jb) = Tsfc0(i,jb)
          Qsfcm(i,jb) = Qsfc0(i,jb)
       END IF
    END DO

    nLndPts=0
    DO i=1,nCols
       IF(iMask(i) >= 1_i8)THEN 
          nLndPts=nLndPts+1
          tgrd(nLndPts)= tg    (nLndPts,jb)
          roff(nLndPts)= grunof(nLndPts,jb)
          ect(nLndPts )= gtrans(nLndPts,jb) * hltm * dtc3x!Transpiracao no topo da copa (J/m*m)
          eci(nLndPts )= ginvap(nLndPts,jb) * hltm * dtc3x!...Evaporacao da agua interceptada no topo da copa (J/m*m) ! (kg m-2 s-1 * J/kg*dt)
          egt(nLndPts )= gtransu(nLndPts,jb)* hltm * dtc3x!Transpiracao na base da copa (J/m*m)
          egi(nLndPts )= gtransl(nLndPts,jb)* hltm * dtc3x !Evaporacao da neve (J/m*m)
          egs(nLndPts )= gsuvap (nLndPts,jb)* hltm * dtc3x !Evaporacao do solo arido (J/m*m)

          !WRITE(*,'(A,5F12.5)')'pkubota', ect(nLndPts ),eci(nLndPts ),egt(nLndPts ),egi(nLndPts ),egs(nLndPts )

          iMaskIBIS(nLndPts,jb) =  INT(vegtype0  (nLndPts,jb)) 
          umom     (i) =  taux  (nLndPts)
          vmom     (i) =  tauy  (nLndPts)
          bstar    (i) =  bstar1(nLndPts)
          sens     (i) =  -fsena (nLndPts)
          evap     (i) =  -fvapa (nLndPts)*hvap
          Tsfc0    (i,jb)       = tg    (nLndPts,jb)
          ! Tsfc0    (i,jb)       = ts2    (nLndPts)
          ! Qsfc0    (i,jb)       = MAX   (1.0e-12_r8,qa (nLndPts))
          Qsfc0    (i,jb)       = MAX   (1.0e-12_r8,qs2 (nLndPts))
          Tsfcm    (i,jb)       = tg    (nLndPts,jb)
          ! Tsfcm    (i,jb)       = ts2    (nLndPts)
          ! Qsfcm    (i,jb)       = MAX   (1.0e-12_r8,qa (nLndPts))       
          Qsfcm    (i,jb)       = MAX   (1.0e-12_r8,qs2 (nLndPts))       
          w0       (nLndPts,1,jb) = wsoi(nLndPts,1,jb)
          w0       (nLndPts,2,jb) = wsoi(nLndPts,3,jb)
          w0       (nLndPts,3,jb) = wsoi(nLndPts,nsoilay,jb)
          wm       (nLndPts,1,jb) = wsoi(nLndPts,1,jb)
          wm       (nLndPts,2,jb) = wsoi(nLndPts,3,jb)
          wm       (nLndPts,3,jb) = wsoi(nLndPts,nsoilay,jb)
          capac0   (nLndPts,1,jb) = wliqu(nLndPts,jb) + wliqs(nLndPts,jb)
          capac0   (nLndPts,2,jb) = wliql(nLndPts,jb)
          capacm   (nLndPts,1,jb) = wliqu(nLndPts,jb) + wliqs(nLndPts,jb)
          capacm   (nLndPts,2,jb) = wliql(nLndPts,jb)
          td0      (nLndPts,jb)   = tsoi(nLndPts,nsoilay,jb)  
          tdm      (nLndPts,jb)   = tsoi(nLndPts,nsoilay,jb)  
          tc0      (nLndPts,jb)   = tu    (nLndPts,jb)
          tcm      (nLndPts,jb)   = tu    (nLndPts,jb)
          tg0      (nLndPts,jb)   = tg    (nLndPts,jb)
          tgm      (nLndPts,jb)   = tg    (nLndPts,jb)
       END IF
    END DO
    !    asurd (1:nLndPts,1:nBand,jb),asuri(1:nLndPts,1:nBand,jb)
    IF(dodia(nDiag_biomau))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             diag(i) = totbiou (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biomau,jb)
    ENDIF
    IF(dodia(nDiag_biomal))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             diag(i) = totbiol (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biomal,jb)
    ENDIF
    IF(dodia(nDiag_tlaiup))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN 
             nLndPts=nLndPts+1
             diag(i) = totlaiu(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tlaiup,jb)
    ENDIF
    IF(dodia(nDiag_tlailw))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = totlail(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tlailw,jb)
    ENDIF
    IF(dodia(nDiag_tstnsp))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = storedn(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tstnsp,jb)
    ENDIF



    IF(dodia(nDiag_somdfa))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = decompl   (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_somdfa,jb)
    ENDIF
    IF(dodia(nDiag_lidecf))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = decomps (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_lidecf,jb)
    ENDIF
    IF(dodia(nDiag_wsttot))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = wtot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_wsttot,jb)
    ENDIF

    IF(dodia(nDiag_facuca))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = fu (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_facuca,jb)
    ENDIF
    IF(dodia(nDiag_fsfclc))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = fl (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_fsfclc,jb)
    ENDIF
    IF(dodia(nDiag_frsnow))THEN !! fractional snow cover
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = fi (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_frsnow,jb)
    ENDIF

    IF(dodia(nDiag_insnpp))THEN ! instantaneous npp (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tnpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_insnpp,jb)
    ENDIF

    IF(dodia(nDiag_insnee))THEN ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tneetot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_insnee,jb)
    ENDIF

    IF(dodia(nDiag_grbdy0))THEN ! annual total growing degree days for current year > 0C
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gdd0this (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_grbdy0,jb)
    ENDIF

    IF(dodia(nDiag_grbdy5))THEN ! annual total growing degree days for current year > 5C  
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gdd5this (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_grbdy5,jb)
    ENDIF

    IF(dodia(nDiag_avet2m))THEN ! monthly average 2-m surface-air temperature 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = amts2 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_avet2m,jb)
    ENDIF

    IF(dodia(nDiag_monnpp))THEN ! monthly total npp for ecosystem (kg-C/m**2/month)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = amnpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_monnpp,jb)
    ENDIF

    IF(dodia(nDiag_monnee))THEN ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = amneetot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_monnee,jb)
    ENDIF

    IF(dodia(nDiag_yeanpp))THEN ! annual total npp for ecosystem (kg-c/m**2/yr)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_yeanpp,jb)
    ENDIF

    IF(dodia(nDiag_yeanee))THEN ! annual total npp for ecosystem (kg-c/m**2/yr)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ayneetot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_yeanee,jb)
    ENDIF




    IF(dodia(nDiag_upclai))THEN ! upper canopy single-sided leaf area index (area leaf/area veg)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = lai (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_upclai,jb)
    ENDIF

    IF(dodia(nDiag_lwclai))THEN ! lower canopy single-sided leaf area index (area leaf/area veg)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = lai (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_lwclai,jb)
    ENDIF

    IF(dodia(nDiag_pfts01))THEN  ! pft tropical broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts01,jb)
    ENDIF

    IF(dodia(nDiag_pfts02))THEN  ! pft tropical broadleaf drought-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts02,jb)
    ENDIF

    IF(dodia(nDiag_pfts03))THEN  ! pft warm-temperate broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,3,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts03,jb)
    ENDIF

    IF(dodia(nDiag_pfts04))THEN  ! pft temperate conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,4,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts04,jb)
    ENDIF

    IF(dodia(nDiag_pfts05))THEN  ! pft temperate broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,5,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts05,jb)
    ENDIF

    IF(dodia(nDiag_pfts06))THEN  ! pft boreal conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,6,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts06,jb)
    ENDIF

    IF(dodia(nDiag_pfts07))THEN  ! pft boreal broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,7,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts07,jb)
    ENDIF

    IF(dodia(nDiag_pfts08))THEN  ! pft boreal conifer cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,8,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts08,jb)
    ENDIF

    IF(dodia(nDiag_pfts09))THEN  ! pft evergreen shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,9,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts09,jb)
    ENDIF

    IF(dodia(nDiag_pfts10))THEN  ! pft cold-deciduous shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,10,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts10,jb)
    ENDIF

    IF(dodia(nDiag_pfts11))THEN  ! pft cool (c4) grasses   
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,11,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts11,jb)
    ENDIF

    IF(dodia(nDiag_pfts12))THEN   ! pft cool (c3) grasses   
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = plai (nLndPts,12,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_pfts12,jb)
    ENDIF

    IF(dodia(nDiag_biol01))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol01,jb)
    ENDIF

    IF(dodia(nDiag_biol02))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol02,jb)
    ENDIF

    IF(dodia(nDiag_biol03))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,3,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol03,jb)
    ENDIF

    IF(dodia(nDiag_biol04))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,4,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol04,jb)
    ENDIF

    IF(dodia(nDiag_biol05))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,5,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol05,jb)
    ENDIF

    IF(dodia(nDiag_biol06))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,6,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol06,jb)
    ENDIF

    IF(dodia(nDiag_biol07))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,7,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol07,jb)
    ENDIF

    IF(dodia(nDiag_biol08))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,8,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol08,jb)
    ENDIF

    IF(dodia(nDiag_biol09))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,9,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol09,jb)
    ENDIF

    IF(dodia(nDiag_biol10))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,10,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol10,jb)
    ENDIF

    IF(dodia(nDiag_biol11))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,11,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol11,jb)
    ENDIF

    IF(dodia(nDiag_biol12))THEN 
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cbiol (nLndPts,12,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_biol12,jb)
    ENDIF

    IF(dodia(nDiag_ynpp01))THEN  ! ynpp tropical broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,1,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp01,jb)
    ENDIF

    IF(dodia(nDiag_ynpp02))THEN ! ynpp tropical broadleaf drought-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,2,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp02,jb)
    ENDIF

    IF(dodia(nDiag_ynpp03))THEN  ! ynpp warm-temperate broadleaf evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,3,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp03,jb)
    ENDIF

    IF(dodia(nDiag_ynpp04))THEN  ! ynpp temperate conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,4,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp04,jb)
    ENDIF

    IF(dodia(nDiag_ynpp05))THEN  ! ynpp temperate broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,5,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp05,jb)
    ENDIF

    IF(dodia(nDiag_ynpp06))THEN  ! ynpp boreal conifer evergreen trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,6,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp06,jb)
    ENDIF

    IF(dodia(nDiag_ynpp07))THEN  !ynpp boreal broadleaf cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,7,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp07,jb)
    ENDIF

    IF(dodia(nDiag_ynpp08))THEN  ! ynpp boreal conifer cold-deciduous trees
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,8,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp08,jb)
    ENDIF

    IF(dodia(nDiag_ynpp09))THEN  ! ynpp evergreen shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,9,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp09,jb)
    ENDIF

    IF(dodia(nDiag_ynpp10))THEN  ! ynpp cold-deciduous shrubs
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,10,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp10,jb)
    ENDIF

    IF(dodia(nDiag_ynpp11))THEN  !! ynpp warm (c4) grasses
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,11,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp11,jb)
    ENDIF

    IF(dodia(nDiag_ynpp12))THEN  ! ynpp cool (c3) grasses
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aynpp (nLndPts,12,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_ynpp12,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cmontp = 170 !coldest monthly temperature                             (C)
    IF(dodia(nDiag_cmontp))THEN  !coldest monthly temperature                             (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cmontp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_wmontp = 171 !warmest monthly temperature                             (C)
    IF(dodia(nDiag_wmontp))THEN   !warmest monthly temperature                             (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tw (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_wmontp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_atogpp = 172 !annual total gpp for ecosystem                               (kg-c/m**2/yr)
    IF(dodia(nDiag_atogpp))THEN  !annual total gpp for ecosystem                               (kg-c/m**2/yr)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = aygpptot (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_atogpp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_toigpp = 173 !instantaneous gpp                                (mol-CO2 / m-2 / second)
    IF(dodia(nDiag_toigpp))THEN  !instantaneous gpp                                (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tgpptot (nLndPts)
          END IF
       END DO
       CALL updia(diag,nDiag_toigpp,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_fxcsol = 174 !instantaneous fine co2 flux from soil                       (mol-CO2 / m-2 / second)
    IF(dodia(nDiag_fxcsol))THEN  !instantaneous fine co2 flux from soil                       (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tco2root(nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_fxcsol,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_mcsoil = 175 !instantaneous microbial co2 flux from soil       (mol-CO2 / m-2 / second)
    IF(dodia(nDiag_mcsoil))THEN  !instantaneous microbial co2 flux from soil       (mol-CO2 / m-2 / second)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tco2mic (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_mcsoil,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcub = 176 !canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcub))THEN !canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcuc = 177 !canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcuc))THEN  !canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcuc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcls = 178 !canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcls))THEN  !canopy average gross photosynthesis rate - shrubs     (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcls,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcl4 = 179 !canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcl4))THEN  !canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcl4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cagcl3 = 180 !canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cagcl3))THEN  !canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = agcl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cagcl3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancub = 181 !canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancub))THEN  !canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancuc = 182 !canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancuc))THEN  !canopy average net photosynthesis rate - conifer      (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancuc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancls = 183 !canopy average net photosynthesis rate - shrubs            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancls))THEN  !canopy average net photosynthesis rate - shrubs            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancls,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancl4 = 184 !canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancl4))THEN  !canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancl4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cancl3 = 185 !canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    IF(dodia(nDiag_cancl3))THEN  !canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ancl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cancl3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicoub = 186 !intercellular co2 concentration - broadleaf            (mol_co2/mol_air)
    IF(dodia(nDiag_cicoub))THEN  !intercellular co2 concentration - broadleaf            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ciub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicoub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicouc = 187 !intercellular co2 concentration - conifer             (mol_co2/mol_air)
    IF(dodia(nDiag_cicouc))THEN  !intercellular co2 concentration - conifer             (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = ciuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicouc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscoub = 188 !leaf boundary layer co2 concentration - broadleaf     (mol_co2/mol_air)
    IF(dodia(nDiag_cscoub))THEN  ! !leaf boundary layer co2 concentration - broadleaf     (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscoub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscoub = 189 !upper canopy stomatal conductance - broadleaf            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscoub))THEN  !upper canopy stomatal conductance - broadleaf            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsub (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscoub,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscouc = 190 !leaf boundary layer co2 concentration - conifer            (mol_co2/mol_air)
    IF(dodia(nDiag_cscouc))THEN   !leaf boundary layer co2 concentration - conifer            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscouc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscouc = 191 !upper canopy stomatal conductance - conifer            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscouc))THEN  !upper canopy stomatal conductance - conifer            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsuc (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscouc,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicols = 192 !intercellular co2 concentration - shrubs              (mol_co2/mol_air)
    IF(dodia(nDiag_cicols))THEN  !intercellular co2 concentration - shrubs              (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cils (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicols,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicol3 = 193 !intercellular co2 concentration - c3 plants            (mol_co2/mol_air)
    IF(dodia(nDiag_cicol3))THEN  !intercellular co2 concentration - c3 plants            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cil3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicol3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cicol4 = 194 !intercellular co2 concentration - c4 plants            (mol_co2/mol_air)
    IF(dodia(nDiag_cicol4))THEN  !intercellular co2 concentration - c4 plants            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = cil4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cicol4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscols = 195 !leaf boundary layer co2 concentration - shrubs            (mol_co2/mol_air)
    IF(dodia(nDiag_cscols))THEN  !leaf boundary layer co2 concentration - shrubs            (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscols,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscols = 196 !lower canopy stomatal conductance - shrubs            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscols))THEN  !lower canopy stomatal conductance - shrubs            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsls (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscols,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscol3 = 197 !leaf boundary layer co2 concentration - c3 plants     (mol_co2/mol_air)
    IF(dodia(nDiag_cscol3))THEN  !leaf boundary layer co2 concentration - c3 plants     (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscol3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscol3 = 198 !lower canopy stomatal conductance - c3 grasses            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscol3))THEN  !lower canopy stomatal conductance - c3 grasses            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsl3 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscol3,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_cscol4 = 199 !leaf boundary layer co2 concentration - c4 plants     (mol_co2/mol_air)
    IF(dodia(nDiag_cscol4))THEN  !leaf boundary layer co2 concentration - c4 plants     (mol_co2/mol_air)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = csl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_cscol4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_gscol4 = 200 !lower canopy stomatal conductance - c4 grasses            (mol_co2 m-2 s-1)
    IF(dodia(nDiag_gscol4))THEN  ! lower canopy stomatal conductance - c4 grasses            (mol_co2 m-2 s-1)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = gsl4 (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_gscol4,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_tcthis = 201 !coldest monthly temperature of current year              (C)
    IF(dodia(nDiag_tcthis))THEN  !coldest monthly temperature of current year              (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = tcthis (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_tcthis,jb)
    ENDIF

    !INTEGER, PUBLIC, PARAMETER :: nDiag_twthis = 201 !warmest monthly temperature of current year              (C)
    IF(dodia(nDiag_twthis))THEN  !warmest monthly temperature of current year              (C)
       diag=0.0_r8
       nLndPts=0
       DO i=1, nCols
          IF(iMask(i) >= 1_i8)THEN  
             nLndPts=nLndPts+1
             diag(i) = twthis (nLndPts,jb)
          END IF
       END DO
       CALL updia(diag,nDiag_twthis,jb)
    ENDIF


    idateprev=idatec
  END SUBROUTINE Ibis_Interface



  SUBROUTINE IbisDrv(mcsec        ,pi          ,stef         ,vonk         ,grav        , &
       tmelt        ,hfus        ,hvap         ,hsub         ,ch2o        , &
       cice         ,cair        ,cvap         ,rair         ,rvap        , &
       cappa        ,rhow        ,npoi         ,nband        ,nsoilay     , & 
       nsnolay      ,npft        ,epsilon      ,dtime        ,doalb       , &
       ginvap       ,gsuvap      ,gtrans       ,gtransu      ,gtransl     , &
       grunof       ,gdrain      ,gadjust      ,a10scalparamu,a10daylightu, &
       a10scalparaml,a10daylightl,vmax_pft     ,tau15        ,kc15        , &
       ko15         ,cimax       ,gammaub      ,alpha3       ,theta3      , &
       beta3        ,coefmub     ,coefbub      ,gsubmin      ,gammauc     , &
       coefmuc      ,coefbuc     ,gsucmin      ,gammals      ,coefmls     , & 
       coefbls      ,gslsmin     ,gammal3      ,coefml3      ,coefbl3     , &
       gsl3min      ,gammal4     ,alpha4       ,theta4       ,beta4       , &
       coefml4      ,coefbl4     ,gsl4min      ,wliqu        ,wliqumax    , & 
       wsnou        ,wsnoumax    ,tu           ,wliqs        ,wliqsmax    , & 
       wsnos        ,wsnosmax    ,ts           ,wliql        ,wliqlmax    , &  
       wsnol        ,wsnolmax    ,tl           ,topparu      ,topparl     , &
       fl           ,fu          ,lai          ,sai          ,rhoveg      , &   
       tauveg       ,orieh       ,oriev        ,wliqmin      ,wsnomin     , & 
       t12          ,tdripu      ,tblowu       ,tdrips       ,tblows      , &
       t34          ,tdripl      ,tblowl       ,ztop         ,alaiml      , &
       zbot         ,alaimu      ,froot        ,q34          ,q12         , &
       su                 ,cleaf       ,dleaf        ,ss           ,cstem       , & 
       dstem        ,sl          ,cgrass       ,ciub         ,ciuc        , &
       exist        ,csub        ,gsub         ,csuc         ,gsuc        , &
       agcub        ,agcuc       ,ancub        ,ancuc        ,totcondub   , &
       totconduc    ,cils        ,cil3         ,cil4         ,csls        , &
       gsls         ,csl3        ,gsl3         ,csl4         ,gsl4        , &
       agcls        ,agcl4       ,agcl3        ,ancls        ,ancl4       , &
       ancl3       ,totcondls    ,totcondl3    ,totcondl4    ,chu         , &
       chs          ,chl         ,frac         ,tlsub          ,z0sno       , & 
       rhos         ,consno      ,hsnotop      ,hsnomin      ,fimin       , &
       fimax        ,fi          ,tsno         ,hsno         ,sand               , &
       clay         ,poros       ,wsoi         ,wisoi        ,consoi      , &  
       zwpmax       ,wpud        ,wipud        ,wpudmax      ,qglif       , &         
       tsoi         ,hvasug      ,hvasui       ,albsav       ,albsan      , &
       tg           ,ti          ,z0soi        ,swilt        ,sfield      , &
       stressl      ,stressu     ,stresstl     ,stresstu     ,csoi        , &         
       rhosoi       ,hsoi        ,suction      ,bex          ,upsoiu      , &  
       upsoil       ,heatg       ,heati        ,hydraul      ,porosflo    , &
       ibex         ,bperm       ,hflo            ,ta           ,asurd       , &  
       asuri        ,coszen      ,solad        ,solai        ,fira               , & 
       raina        ,qa          ,psurf        ,snowa        ,ua               , &   
       o2conc       ,co2conc     ,&
       td                 ,vzero       ,ndaypy       ,nppdummy     , & 
       cbiow        ,sapfrac      ,cbior       , & 
       tco2root    ,tneetot      ,tco2mic      ,a10td       , &
       a10ancub     ,a10ancuc    ,a10ancls     ,a10ancl3     ,a10ancl4    , & 
       ndtimes      ,adrain       ,adsnow      ,tnpptot, &
       adaet        ,adtrunoff   ,adsrunoff    ,addrainage   ,adrh        , &
       adsnod       ,adsnof      ,adwsoi       ,adtsoi       ,adwisoi     , &
       adtlaysoi    ,adwlaysoi   ,adwsoic      ,adtsoic      ,adco2mic    , &
       adco2root    ,adco2soi    ,adco2ratio   ,adnmintot    ,decompl     , &    
       decomps      ,tnmin       ,ndaypm       ,nmtimes     , & 
       amrain       ,amsnow      ,amaet        ,amtrunoff    ,amsrunoff   , &
       amdrainage   ,amtemp      ,amqa         , &
       amsolar      ,amirup      ,amirdown     ,amsens       ,amlatent    , &  
       amlaiu       ,amlail      ,amtsoi       ,amwsoi       ,amwisoi     , &   
       amvwc        ,amawc       ,amsnod       ,amsnof       ,amnpp       , &
       amnpptot     ,amco2mic    ,amco2root    ,amco2soi     ,amco2ratio  , &
       amneetot     ,amnmintot   ,amalbedo     ,amtsoil      ,amwsoil     , & 
       amwisoil     ,nytimes      ,aysolar      ,ayirup      , &
       ayirdown     ,aysens      ,aylatent     ,ayprcp       ,ayaet       , &  
       aytrans      ,aytrunoff   ,aysrunoff    ,aydrainage   ,aydwtot     , & 
       aywsoi       ,aywisoi     ,aytsoi       ,ayvwc        ,ayawc       , &  
       aystresstu   ,aystresstl  ,aygpp            ,aygpptot     ,aynpp       , & 
       aynpptot     ,ayco2mic    ,ayco2root    ,ayco2soi     ,ayneetot    , &
       ayrootbio    ,aynmintot   ,ayalit       ,ayblit       ,aycsoi      , & 
       aycmic       ,ayanlit     ,aybnlit      ,aynsoi       ,ayalbedo    ,&
       totalit     , &
       totrlit      ,totcsoi     ,totcmic      ,totanlit     ,totrnlit    , &
       totnsoi      ,totnmic     ,totlit       ,totfall      ,totnlit     , &
       firefac      ,wtot        ,storedn      ,yrleach      ,ynleach     , & 
       falll        ,fallr       ,fallw        ,clitlm       ,clitls      , &
       clitrm       ,clitrs      ,clitwm       ,clitws       ,csoislop    , &
       csoislon     ,csoipas     ,clitll       ,clitrl       ,clitwl      , &  
       tc           ,agddu              ,tempu            ,agddl        ,templ       , &
       dropu        ,dropls      ,dropl4       ,dropl3       ,plai               , &
       iday        ,imonth       ,iyear        ,iyear0       , &
       isimveg      ,spinmax     ,amts2        ,amtransu     ,amtransl    , &
       amsuvap      ,aminvap     ,ux            ,uy                  ,taux        , &
       tauy         ,ts2              ,qs2          ,deltat       ,gdd0        , &  
       gdd0this     ,tcthis      ,twthis       ,tcmin        ,gdd5        , & 
       gdd5this     ,TminL       ,TminU        ,Twarm        ,GDD         , & 
       aleaf        ,awood       ,cbiol        ,aroot        ,disturbf    , &
       disturbo     ,specla      ,biomass      ,totlaiu      ,totlail     , &  
       totbiou      ,totbiol     ,woodnorm     ,vegtype0     ,tauwood0    , & 
       tauwood      ,tauleaf     ,tauroot      ,xminlai      ,cdisturb    , & 
       ayanpp       ,ayanpptot   , &
       nstep         ,imonthprev   , &
       isimco2      ,isimfire    , &
       co2init      ,tw           ,fvapa       ,fsena        , &
       z0           ,ustar       ,hc           ,hg          ,ec           , &
       eg           ,dispu       ,cu           ,firb        ,tgpptot      , &
       bstar)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN        ) :: mcsec!global  ! current seconds in day (0 - (86400 - dtime))
    REAL(KIND=r8), INTENT(IN        ) :: pi   !global
    REAL(KIND=r8), INTENT(IN        ) :: stef !global  ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8), INTENT(IN        ) :: vonk !global  ! von karman constant (dimensionless)
    REAL(KIND=r8), INTENT(IN        ) :: grav !global  ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN        ) :: tmelt!global  ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN        ) :: hfus !global  ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hvap !global  ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hsub !global  ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: ch2o !global  ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cice !global  ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cair !global  ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cvap !global  ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rair !global  ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rvap !global  ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cappa!global  ! rair/cair
    REAL(KIND=r8), INTENT(IN        ) :: rhow !global  ! density of liquid water (all types) (kg m-3)

    ! 
    !
    INTEGER, INTENT(IN   ) :: npoi   !global  
    INTEGER, INTENT(IN   ) :: nband  !global  
    INTEGER, INTENT(IN   ) :: nsoilay!global   ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay!global   ! number of snow layers
    INTEGER, INTENT(IN   ) :: npft   !global   ! number of plant functional types
    REAL(KIND=r8), INTENT(IN   ) :: epsilon!global   ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   )  :: dtime !global   ! model timestep (seconds)
    LOGICAL, INTENT(IN   )  :: doalb !global    ! true if surface albedo calculation time step

    !      INCLUDE 'comhyd.h'
    REAL(KIND=r8), INTENT(OUT  ) :: ginvap (npoi)!local ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gsuvap (npoi)!local ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gtrans (npoi)!local ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransu(npoi)!local ! transpiration from upper canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransl(npoi)!local ! transpiration from lower canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: grunof (npoi)!local ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gdrain (npoi)!local ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gadjust(npoi)!local ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
    !      INCLUDE 'comsum.h'
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparamu(npoi)!global ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightu (npoi)!global ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparaml(npoi)!global ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightl (npoi)!global ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
    !      INCLUDE 'compft.h'
    REAL(KIND=r8), INTENT(IN   ) :: vmax_pft(npft)!global ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]
    REAL(KIND=r8), INTENT(IN   ) :: tau15           !global ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: kc15           !global ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: ko15           !global ! o2 kinetic parameter (mol/mol) 
    REAL(KIND=r8), INTENT(IN   ) :: cimax           !global ! maximum value for ci (needed for model stability)
    REAL(KIND=r8), INTENT(IN   ) :: gammaub           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha3           !global ! intrinsic quantum efficiency for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: beta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefmub           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbub           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsubmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammauc           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmuc           !global ! 'm' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: coefbuc           !global ! 'b' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: gsucmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammals           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmls           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbls           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gslsmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal3           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefml3           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbl3           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gsl3min           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal4           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha4           !global ! intrinsic quantum efficiency for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless) 
    REAL(KIND=r8), INTENT(IN   ) :: beta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefml4           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbl4           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsl4min           !global ! absolute minimum stomatal conductance
    !      include 'comveg.h'
    REAL(KIND=r8), INTENT(INOUT) :: wliqu    (npoi)  !global ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqumax              !global ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnou    (npoi)  !global ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnoumax              !global ! intercepted snow capacity for upper canopy leaves (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tu       (npoi)  !global ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqs    (npoi)  !global ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqsmax              !global ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnos    (npoi)  !global ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnosmax              !global ! intercepted snow capacity for upper canopy stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: ts       (npoi)  !global ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliql    (npoi)  !global ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqlmax              !global ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnol    (npoi)  !global ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnolmax              !global ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tl       (npoi)  !global ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(INOUT) :: topparu  (npoi)  !local  ! total photosynthetically active raditaion absorbed 
    ! by top leaves of upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: topparl  (npoi)  !local  ! total photosynthetically active raditaion absorbed
    ! by top leaves of lower canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: fl       (npoi)   !global ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(INOUT) :: fu       (npoi)   !global ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: lai      (npoi,2) !global ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(INOUT) :: sai      (npoi,2) !global ! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: rhoveg   (nband,2)!global ! reflectance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: tauveg   (nband,2)!global  ! transmittance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: orieh    (2)      !global! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8), INTENT(IN   ) :: oriev    (2)      !global! fraction of leaf/stems with vertical
    REAL(KIND=r8), INTENT(INOUT) :: wliqmin           !local ! minimum intercepted water on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnomin           !local ! minimum intercepted snow on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: t12      (npoi)   !global ! air temperature at z12 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripu            !global ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowu             !global ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tdrips             !global ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
    REAL(KIND=r8), INTENT(IN   ) :: tblows             !global ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: t34      (npoi)   !global ! air temperature at z34 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripl            !global ! decay time for dripoff of liquid intercepted
    ! by lower canopy leaves & stem (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowl             ! global          ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: ztop     (npoi,2) ! global  ! height of plant top above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaiml             ! global ! lower canopy leaf & stem maximum area (2 sided) for
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: zbot     (npoi,2) ! global  ! height of lowest branches above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaimu             ! global  ! upper canopy leaf & stem area (2 sided) for 
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: froot    (nsoilay,2)! global! fraction of root in soil layer 
    REAL(KIND=r8), INTENT(INOUT) :: q34      (npoi)   ! global! specific humidity of air at z34
    REAL(KIND=r8), INTENT(INOUT) :: q12      (npoi)   ! global! specific humidity of air at z12
    REAL(KIND=r8), INTENT(INOUT) :: su       (npoi)   ! local ! air-vegetation transfer coefficients (*rhoa) for
    ! upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cleaf             ! global! empirical constant in upper canopy leaf-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dleaf    (2)             ! global ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: ss       (npoi)   ! local! air-vegetation transfer coefficients (*rhoa) for 
    ! upper canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cstem             ! global ! empirical constant in upper canopy stem-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dstem    (2)             ! global ! typical linear stem dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: sl       (npoi)   ! local ! air-vegetation transfer coefficients (*rhoa) for 
    ! lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cgrass             ! global ! empirical constant in lower canopy-air aerodynamic 
    ! transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(INOUT) :: ciub     (npoi)         ! global ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: ciuc     (npoi)         ! global ! intercellular co2 concentration - conifer        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: exist    (npoi,npft)  ! global ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(INOUT) :: csub     (npoi)         ! global ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsub     (npoi)         ! global ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csuc     (npoi)         ! global ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsuc     (npoi)         ! global ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcub    (npoi)         ! local  ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcuc    (npoi)         ! local  ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancub    (npoi)         ! local  ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancuc    (npoi)         ! local  ! canopy average net photosynthesis rate - conifer          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondub(npoi)         ! local  ! 
    REAL(KIND=r8), INTENT(INOUT) :: totconduc(npoi)         ! local  !
    REAL(KIND=r8), INTENT(INOUT) :: cils     (npoi)         ! global ! intercellular co2 concentration - shrubs        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil3     (npoi)         ! global ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil4     (npoi)         ! global ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: csls     (npoi)         ! global ! leaf boundary layer co2 concentration - shrubs   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsls     (npoi)         ! global ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl3     (npoi)         ! global ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl3     (npoi)         ! global ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl4     (npoi)         ! global ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl4     (npoi)         ! global ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcls    (npoi)         ! local  ! canopy average gross photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl4    (npoi)         ! local ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl3    (npoi)         ! local ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancls    (npoi)         ! local ! canopy average net photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl4    (npoi)         ! local ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl3    (npoi)         ! local ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondls(npoi)         ! local ! 
    REAL(KIND=r8), INTENT(INOUT) :: totcondl3(npoi)         ! local !
    REAL(KIND=r8), INTENT(INOUT) :: totcondl4(npoi)         ! local !
    REAL(KIND=r8), INTENT(IN   ) :: chu                         ! global ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chs                         ! global ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chl                         ! global ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: frac     (npoi,npft)  ! global ! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8), INTENT(INOUT) :: tlsub    (npoi)         ! global ! temperature of lower canopy vegetation buried by snow (K)
    !      INCLUDE 'comsat.h'    
    !      include 'comsno.h'
    REAL(KIND=r8), INTENT(IN   ) :: z0sno  ! global ! roughness length of snow surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: rhos   ! global ! density of snow (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: consno ! global ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: hsnotop! global ! thickness of top snow layer (m)
    REAL(KIND=r8), INTENT(IN   ) :: hsnomin! global ! minimum total thickness of snow (m)
    REAL(KIND=r8), INTENT(IN   ) :: fimin  ! global ! minimum fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: fimax  ! global ! maximum fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: fi     (npoi)! global ! fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: tsno   (npoi,nsnolay)! global ! temperature of snow layers (K)
    REAL(KIND=r8), INTENT(INOUT) :: hsno   (npoi,nsnolay)! global ! thickness of snow layers (m)

    !      INCLUDE 'comsoi.h'
    REAL(KIND=r8), INTENT(IN   ) :: sand    (npoi,nsoilay)! global ! percent sand of soil
    REAL(KIND=r8), INTENT(IN   ) :: clay    (npoi,nsoilay)! global ! percent clay of soil
    REAL(KIND=r8), INTENT(IN   ) :: poros   (npoi,nsoilay)! global ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8), INTENT(INOUT) :: wsoi    (npoi,nsoilay)! global ! fraction of soil pore space containing liquid water
    REAL(KIND=r8), INTENT(INOUT) :: wisoi   (npoi,nsoilay)! global ! fraction of soil pore space containing ice
    REAL(KIND=r8), INTENT(INOUT) :: consoi  (npoi,nsoilay)! local  ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: zwpmax                 ! global ! assumed maximum fraction of soil surface 
    ! covered by puddles (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: wpud    (npoi)! global ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wipud   (npoi)! global ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wpudmax         ! global ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: qglif   (npoi,4) ! local ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    REAL(KIND=r8), INTENT(INOUT) :: tsoi    (npoi,nsoilay)! global        ! soil temperature for each layer (K)
    REAL(KIND=r8), INTENT(INOUT) :: hvasug  (npoi)        ! local ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8), INTENT(INOUT) :: hvasui  (npoi)        ! local ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: albsav  (npoi)        ! global ! saturated soil surface albedo (visible waveband)
    REAL(KIND=r8), INTENT(IN   ) :: albsan  (npoi)        ! global ! saturated soil surface albedo (near-ir waveband)
    REAL(KIND=r8), INTENT(INOUT) :: tg      (npoi)        ! global ! soil skin temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: ti      (npoi)        ! global ! snow skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: z0soi   (npoi)        ! global ! roughness length of soil surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: swilt   (npoi,nsoilay)! global ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(IN   ) :: sfield  (npoi,nsoilay)! global ! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(INOUT) :: stressl (npoi,nsoilay)! local ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stressu (npoi,nsoilay)! local ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstl(npoi)        ! local ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstu(npoi)        ! local ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: csoi    (npoi,nsoilay)! global ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8), INTENT(IN   ) :: rhosoi  (npoi,nsoilay)! global ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: hsoi    (nsoilay+1)   ! global ! soil layer thickness (m)
    REAL(KIND=r8), INTENT(IN   ) :: suction (npoi,nsoilay)! global ! saturated matric potential (m-h2o)
    REAL(KIND=r8), INTENT(IN   ) :: bex     (npoi,nsoilay)! global ! exponent "b" in soil water potential
    REAL(KIND=r8), INTENT(INOUT) :: upsoiu  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: upsoil  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: heatg   (npoi)         ! local  ! net heat flux into soil surface (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: heati   (npoi)         ! local  ! net heat flux into snow surface (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: hydraul (npoi,nsoilay)! global ! saturated hydraulic conductivity (m/s)
    REAL(KIND=r8), INTENT(INOUT) :: porosflo(npoi,nsoilay)! global ! porosity after reduction by ice content
    INTEGER      , INTENT(IN   ) :: ibex    (npoi,nsoilay)! global ! nint(bex), used for cpu speed
    REAL(KIND=r8), INTENT(IN   ) :: bperm                 ! global ! lower b.c. for soil profile drainage 
    ! (0.0 = impermeable; 1.0 = fully permeable)
    REAL(KIND=r8), INTENT(INOUT) :: hflo    (npoi,nsoilay+1)  ! downward heat transport through soil layers (W m-2)


    !   INCLUDE 'comatm.h'
    REAL(KIND=r8), INTENT(IN   ) :: ta     (npoi)         ! global ! air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: asurd  (npoi,nband)   ! local  ! direct albedo of surface system
    REAL(KIND=r8), INTENT(INOUT) :: asuri  (npoi,nband)   ! local  ! diffuse albedo of surface system 
    REAL(KIND=r8), INTENT(IN   ) :: coszen (npoi)         ! global ! cosine of solar zenith angle
    REAL(KIND=r8), INTENT(IN   ) :: solad  (npoi,nband)   ! global ! direct downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: solai  (npoi,nband)   ! global ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fira   (npoi)         ! global ! incoming ir flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: raina  (npoi)         ! global ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: qa     (npoi)         ! global ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8), INTENT(IN   ) :: psurf  (npoi)         ! global ! surface pressure (Pa)
    REAL(KIND=r8), INTENT(IN   ) :: snowa  (npoi)         ! global ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8), INTENT(IN   ) :: ua     (npoi)         ! global ! wind speed (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: o2conc                 ! global ! o2 concentration (mol/mol)
    REAL(KIND=r8), INTENT(INOUT) :: co2conc                 ! global ! co2 concentration (mol/mol)
    REAL(KIND=r8), INTENT(INOUT) :: z0(npoi)
    REAL(KIND=r8), INTENT(INOUT) :: ustar(npoi)
    REAL(KIND=r8), INTENT(OUT) :: hc (npoi)
    REAL(KIND=r8), INTENT(OUT) :: hg (npoi)
    REAL(KIND=r8), INTENT(OUT) :: ec (npoi)
    REAL(KIND=r8), INTENT(OUT) :: eg (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: dispu    (npoi)          ! local ! zero-plane displacement height for upper canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: cu       (npoi)          ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for
    !         upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: firb     (npoi)          ! local ! net upward ir radiation at reference
    ! atmospheric level za (W m-2)
    REAL(KIND=r8), INTENT(OUT) :: bstar(npoi)

    !   INCLUDE 'com1d.h'
    REAL(KIND=r8) :: fwetu    (npoi)          ! local ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqu    (npoi)          ! local ! proportion of fwetu due to liquid
    REAL(KIND=r8) :: fwets    (npoi)          ! local ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqs    (npoi)          ! local ! proportion of fwets due to liquid
    REAL(KIND=r8) :: fwetl    (npoi)          ! local ! fraction of lower canopy stem & leaf area wetted by
    !         intercepted liquid and/or snow
    REAL(KIND=r8) :: rliql    (npoi)          ! local ! proportion of fwetl due to liquid
    REAL(KIND=r8) :: solu     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by upper 
    !         canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8) :: sols     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by upper 
    !         canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: soll     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by lower 
    !         canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: solg     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by unit 
    !         snow-free soil (W m-2)
    REAL(KIND=r8) :: soli     (npoi)          ! local ! solar flux (direct + diffuse) absorbed by unit 
    ! snow surface (W m-2)
    REAL(KIND=r8) :: scalcoefl(npoi,4)     ! local ! term needed in lower canopy scaling
    REAL(KIND=r8) :: scalcoefu(npoi,4)     ! local ! term needed in upper canopy scaling
    INTEGER :: indsol   (npoi)                  ! local ! index of current strip for points with positive coszen
    REAL(KIND=r8) :: albsod   (npoi)          ! local ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsoi   (npoi)          ! local ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsnd   (npoi)          ! local ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: albsni   (npoi)          ! local ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: relod    (npoi)          ! local ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
    REAL(KIND=r8) :: reloi    (npoi)          ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: reupd    (npoi)          ! local ! upward direct radiation per unit incident direct 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: reupi    (npoi)          ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: ablod    (npoi)          ! local ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8) :: abloi    (npoi)          ! local ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8) :: flodd    (npoi)          ! local ! downward direct radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: dummy    (npoi)          ! local ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8) :: flodi    (npoi)          ! local ! downward diffuse radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: floii    (npoi)          ! local ! downward diffuse radiation per unit incident 
    ! diffuse radiation on lower canopy
    REAL(KIND=r8) :: terml    (npoi,7)     ! local ! term needed in lower canopy scaling
    REAL(KIND=r8) :: termu    (npoi,7)     ! local ! term needed in upper canopy scaling
    REAL(KIND=r8) :: abupd    (npoi)          ! local ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8) :: abupi    (npoi)          ! local ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8) :: fupdd    (npoi)          ! local ! downward direct radiation per unit incident direct
    ! beam on upper canopy (W m-2)
    REAL(KIND=r8) :: fupdi    (npoi)          ! local ! downward diffuse radiation per unit icident direct
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: fupii    (npoi)          ! local ! downward diffuse radiation per unit incident diffuse
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: sol2d    (npoi)          ! local ! direct downward radiation  out of upper canopy 
    ! per unit vegetated (upper) area (W m-2)
    REAL(KIND=r8) :: sol2i    (npoi)          ! local ! diffuse downward radiation out of upper
    ! canopy per unit vegetated (upper) area(W m-2)
    REAL(KIND=r8) :: sol3d    (npoi)          ! local ! direct downward radiation  out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8) :: sol3i    (npoi)          ! local ! diffuse downward radiation out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8) :: firs     (npoi)          ! local ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8) :: firu     (npoi)          ! local ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8) :: firl     (npoi)          ! local ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8) :: firg     (npoi)          ! local ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8) :: firi     (npoi)          ! local ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8) :: snowg    (npoi)          ! local ! snowfall rate at soil level (kg h2o m-2 s-1)
    REAL(KIND=r8) :: tsnowg   (npoi)          ! local ! snowfall temperature at soil level (K) 
    REAL(KIND=r8) :: tsnowl   (npoi)          ! local ! snowfall temperature below upper canopy (K)
    REAL(KIND=r8) :: pfluxl   (npoi)          ! local ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) :: raing    (npoi)          ! local ! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8) :: traing   (npoi)          ! local ! rainfall temperature at soil level (K)
    REAL(KIND=r8) :: trainl   (npoi)          ! local ! rainfall temperature below upper canopy (K)
    REAL(KIND=r8) :: snowl    (npoi)          ! local ! snowfall rate below upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8) :: tsnowu   (npoi)          ! local ! snowfall temperature above upper canopy (K)
    REAL(KIND=r8) :: pfluxu   (npoi)          ! local ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8) :: rainu    (npoi)          ! local ! rainfall rate above upper canopy (kg m-2 s-1)
    REAL(KIND=r8) :: trainu   (npoi)          ! local ! rainfall temperature above upper canopy (K)
    REAL(KIND=r8) :: snowu    (npoi)          ! local ! snowfall rate above upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8) :: pfluxs   (npoi)          ! local ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8) :: rainl    (npoi)          ! local ! rainfall rate below upper canopy (kg m-2 s-1)
    REAL(KIND=r8) :: tfac                  ! local ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8) :: cp       (npoi)          ! local ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8) :: za       (npoi)          ! local ! height above the surface of atmospheric forcing (m)
    REAL(KIND=r8) :: bdl      (npoi)          ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for
    ! laower canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) :: dil      (npoi)          ! local ! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8) :: z3       (npoi)          ! local ! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8) :: z4       (npoi)          ! local ! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8) :: z34      (npoi)          ! local ! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8) :: exphl    (npoi)          ! local ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) :: expl     (npoi)          ! local ! exphl**2
    REAL(KIND=r8) :: displ    (npoi)          ! local ! zero-plane displacement height for lower canopy (m)
    REAL(KIND=r8) :: bdu      (npoi)          ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for upper
    ! canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8) :: diu      (npoi)          ! local ! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8) :: z1       (npoi)          ! local ! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8) :: z2       (npoi)          ! local ! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8) :: z12      (npoi)          ! local ! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8) :: exphu    (npoi)          ! local ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8) :: expu     (npoi)          ! local ! exphu**2
    REAL(KIND=r8) :: alogg    (npoi)          ! local ! log of soil roughness
    REAL(KIND=r8) :: alogi    (npoi)          ! local ! log of snow roughness
    REAL(KIND=r8) :: alogav   (npoi)          ! local ! average of alogi and alogg 
    REAL(KIND=r8) :: alog4    (npoi)          ! local ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8) :: alog3    (npoi)          ! local ! log (z3 - displ)
    REAL(KIND=r8) :: alog2    (npoi)          ! local ! log (z2 - displ)
    REAL(KIND=r8) :: alog1    (npoi)          ! local ! log (z1 - dispu) 
    REAL(KIND=r8) :: aloga    (npoi)          ! local ! log (za - dispu) 
    REAL(KIND=r8) :: u2       (npoi)          ! local ! wind speed at level z2 (m s-1)
    REAL(KIND=r8) :: alogu    (npoi)          ! local ! log (roughness length of upper canopy)
    REAL(KIND=r8) :: alogl    (npoi)          ! local ! log (roughness length of lower canopy)
    REAL(KIND=r8) :: richl    (npoi)          ! local ! richardson number for air above upper canopy (z3 to z2)
    REAL(KIND=r8) :: straml   (npoi)          ! local ! momentum correction factor for stratif between
    ! upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8) :: strahl   (npoi)          ! local ! heat/vap correction factor for stratif between
    !         upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8) :: richu    (npoi)          ! local ! richardson number for air between upper & lower canopy (z1 to za)
    REAL(KIND=r8) :: stramu   (npoi)          ! local ! momentum correction factor for stratif above
    !         upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8) :: strahu   (npoi)          ! local ! heat/vap correction factor for stratif above
    !         upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8) :: u1       (npoi)          ! local ! wind speed at level z1 (m s-1)
    REAL(KIND=r8) :: u12      (npoi)          ! local ! wind speed at level z12 (m s-1)
    REAL(KIND=r8) :: u3       (npoi)          ! local ! wind speed at level z3 (m s-1)
    REAL(KIND=r8) :: u34      (npoi)          ! local ! wind speed at level z34 (m s-1)
    REAL(KIND=r8) :: u4       (npoi)          ! local ! wind speed at level z4 (m s-1)
    REAL(KIND=r8) :: cl       (npoi)          ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3)
    !         between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8) :: sg       (npoi)          ! local ! air-soil transfer coefficient
    REAL(KIND=r8) :: si       (npoi)          ! local ! air-snow transfer coefficient
    REAL(KIND=r8) :: fwetux   (npoi)          ! local ! fraction of upper canopy leaf area wetted if dew forms
    REAL(KIND=r8) :: fwetsx   (npoi)          ! local ! fraction of upper canopy stem area wetted if dew forms
    REAL(KIND=r8) :: fwetlx   (npoi)          ! local ! fraction of lower canopy leaf and stem area wetted if dew forms
    REAL(KIND=r8) :: fsena    (npoi)          ! local ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8) :: fseng    (npoi)          ! local ! upward sensible heat flux between soil surface & air at z34 (W m-2)
    REAL(KIND=r8) :: fseni    (npoi)          ! local ! upward sensible heat flux between snow surface & air at z34 (W m-2)
    REAL(KIND=r8) :: fsenu    (npoi)          ! local ! sensible heat flux from upper canopy leaves to air (W m-2)
    REAL(KIND=r8) :: fsens    (npoi)          ! local ! sensible heat flux from upper canopy stems to air (W m-2)
    REAL(KIND=r8) :: fsenl    (npoi)          ! local ! sensible heat flux from lower canopy to air (W m-2)
    REAL(KIND=r8) :: fvapa    (npoi)          ! local ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8) :: fvaput   (npoi)          ! local ! h2o vapor flux (transpiration from dry parts) 
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8) :: fvaps    (npoi)       ! local ! h2o vapor flux (evaporation from wet surface)
    !            between upper canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8) :: fvaplw   (npoi)       ! local ! h2o vapor flux (evaporation from wet surface) 
    !            between lower canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8) :: fvaplt   (npoi)       ! local ! h2o vapor flux (transpiration) 
    !            between lower canopy & air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
    REAL(KIND=r8) :: fvapg    (npoi)       ! local ! h2o vapor flux (evaporation) between soil & air 
    !         at z34 (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8) :: fvapi    (npoi)       ! local ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8) :: fvapuw   (npoi)       ! local ! h2o vapor flux (evaporation from wet parts)
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(INOUT) :: td (npoi)   ! global! daily average temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: vzero(npoi) ! global! a real array of zeros, of length npoi

    INTEGER, INTENT(IN   ) :: ndaypy               ! global! number of days per year
    REAL(KIND=r8), INTENT(OUT  ) :: nppdummy (npoi,npft)! local ! canopy NPP before accounting for stem and root respiration
    REAL(KIND=r8) :: tgpp     (npoi,npft)                 ! local ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tgpptot  (npoi)                         ! local ! instantaneous gpp (mol-CO2 / m-2 / second)
    REAL(KIND=r8) :: tnpp     (npoi,npft)                 ! local ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: cbiow    (npoi,npft)  ! global! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: sapfrac  (npoi)         ! global! fraction of woody biomass that is in sapwood
    REAL(KIND=r8), INTENT(INOUT) :: cbior    (npoi,npft)  ! global! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: tnpptot  (npoi)                         ! local ! instantaneous npp (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tco2root (npoi)         ! local ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tneetot  (npoi)         ! local ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
    REAL(KIND=r8), INTENT(INOUT) :: tco2mic  (npoi)         ! local ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: a10td    (npoi)       ! global! 10-day average daily air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancub (npoi)       ! global! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancuc (npoi)       ! global! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancls (npoi)       ! global! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl3 (npoi)       ! global! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl4 (npoi)       ! global! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)

    INTEGER, INTENT(INOUT) :: ndtimes        (npoi)! global! counter for daily average calculations
    REAL(KIND=r8), INTENT(INOUT) :: adrain    (npoi)! global! daily average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsnow    (npoi)! global! daily average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adaet     (npoi)! global! daily average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adtrunoff (npoi)! global! daily average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsrunoff (npoi)! global! daily average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: addrainage(npoi)! global! daily average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adrh      (npoi)! global! daily average rh (percent)
    REAL(KIND=r8), INTENT(INOUT) :: adsnod    (npoi)! global! daily average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: adsnof    (npoi)! global! daily average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoi    (npoi)! global! daily average soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoi    (npoi)! global! daily average soil temperature (c)
    REAL(KIND=r8), INTENT(INOUT) :: adwisoi   (npoi)! global! daily average soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtlaysoi (npoi)! global! daily average soil temperature (c) of top layer
    REAL(KIND=r8), INTENT(INOUT) :: adwlaysoi (npoi)! global! daily average soil moisture of top layer(fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoic   (npoi)! global! daily average soil moisture using root profile weighting (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoic   (npoi)! global! daily average soil temperature (c) using profile weighting
    REAL(KIND=r8), INTENT(INOUT) :: adco2mic  (npoi)! global! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2root (npoi)! global! daily accumulated co2 respiration from roots (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2soi  (npoi)! global! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2ratio(npoi)! global! ratio of root to total co2 respiration
    REAL(KIND=r8), INTENT(INOUT) :: adnmintot (npoi)! global! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: decompl   (npoi)! global! litter decomposition factor              (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: decomps   (npoi)! global! soil organic matter decomposition factor      (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: tnmin     (npoi)! global! instantaneous nitrogen mineralization (kg_N m-2/timestep)

    INTEGER, INTENT(IN   ) :: ndaypm          (12)          ! global! number of days per month

    INTEGER, INTENT(INOUT) :: nmtimes        (npoi)                   ! global! counter for monthly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: amrain        (npoi)     ! global! monthly average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsnow        (npoi)     ! global! monthly average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amaet        (npoi)     ! global! monthly average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtrunoff    (npoi)     ! global! monthly average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsrunoff    (npoi)     ! global! monthly average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amdrainage   (npoi)     ! global! monthly average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtemp        (npoi)     ! global! monthly average air temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amqa         (npoi)     ! global! monthly average specific humidity (kg-h2o/kg-air)
    REAL(KIND=r8), INTENT(INOUT) :: amsolar        (npoi)     ! global! monthly average incident solar radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirup        (npoi)     ! global! monthly average upward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirdown     (npoi)     ! global! monthly average downward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amsens        (npoi)     ! global! monthly average sensible heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlatent     (npoi)     ! global! monthly average latent heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlaiu        (npoi)     ! global! monthly average lai for upper canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlail        (npoi)     ! global! monthly average lai for lower canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amtsoi        (npoi)     ! global! monthly average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amwsoi        (npoi)     ! global! monthly average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amwisoi        (npoi)     ! global! monthly average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amvwc        (npoi)     ! global! monthly average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amawc        (npoi)     ! global! monthly average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amsnod        (npoi)     ! global! monthly average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: amsnof        (npoi)     ! global! monthly average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amnpp        (npoi,npft)! global! monthly total npp for each plant type (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amnpptot     (npoi)     ! local ! monthly total npp for ecosystem (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2mic     (npoi)     ! global! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2root    (npoi)     ! global! monthly total CO2 flux from soil due to root
    ! respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2soi     (npoi)     ! local ! monthly total soil CO2 flux from microbial
    ! and root respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2ratio   (npoi)       ! local ! monthly ratio of root to total co2 flux
    REAL(KIND=r8), INTENT(OUT  ) :: amneetot     (npoi)       ! local ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amnmintot    (npoi)       ! global! monthly total N mineralization from microbes (kg-N/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amts2        (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransu     (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransl     (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amsuvap      (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: aminvap      (npoi)       ! global
    REAL(KIND=r8), INTENT(INOUT) :: amalbedo     (npoi)         
    REAL(KIND=r8), INTENT(INOUT) :: amtsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwisoil   (npoi, nsoilay)

    INTEGER, INTENT(INOUT) :: nytimes (npoi)                ! global! counter for yearly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: aysolar        (npoi)     ! global! annual average incident solar radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirup        (npoi)     ! global! annual average upward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirdown   (npoi)       ! global! annual average downward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aysens        (npoi)     ! global! annual average sensible heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aylatent   (npoi)       ! global! annual average latent heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayprcp        (npoi)     ! global! annual average precipitation (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayaet        (npoi)     ! global! annual average aet (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrans        (npoi)     ! global! annual average transpiration (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrunoff  (npoi)       ! global! annual average total runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aysrunoff  (npoi)       ! global! annual average surface runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydrainage (npoi)       ! global! annual average drainage (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydwtot        (npoi)     ! global! annual average soil+vegetation+snow water 
    ! recharge (mm/yr or kg_h2o/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aywsoi        (npoi)     ! global! annual average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aywisoi        (npoi)     ! global! annual average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aytsoi        (npoi)     ! global! annual average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: ayvwc        (npoi)     ! global! annual average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: ayawc        (npoi)     ! global! annual average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstu (npoi)       ! global! annual average soil moisture stress 
    ! parameter for upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstl(npoi)        ! global! annual average soil moisture stress 
    ! parameter for lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aygpp     (npoi,npft)   ! global! annual gross npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aygpptot  (npoi)        ! local ! annual total gpp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aynpp     (npoi,npft)   ! global! annual total npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aynpptot  (npoi)        ! local ! annual total npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2mic  (npoi)        ! global! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2root (npoi)        ! global! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayco2soi  (npoi)        ! local ! annual total soil CO2 flux from microbial and 
    ! root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayneetot  (npoi)        ! local ! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayrootbio (npoi)        ! global! annual average live root biomass (kg-C / m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynmintot (npoi)        ! global! annual total nitrogen mineralization (kg-N/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayalit    (npoi)        ! global! aboveground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayblit    (npoi)        ! global! belowground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycsoi    (npoi)        ! global! total soil carbon (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycmic    (npoi)        ! global! total soil carbon in microbial biomass (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayanlit   (npoi)        ! global! aboveground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aybnlit   (npoi)        ! global! belowground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynsoi    (npoi)        ! global! total soil nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayalbedo  (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: totalit   (npoi)           ! global! total standing aboveground litter (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrlit   (npoi)           ! global! total root litter carbon belowground (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcsoi   (npoi)           ! global! total carbon in all soil pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcmic   (npoi)           ! global! total carbon residing in microbial pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totanlit  (npoi)           ! global! total standing aboveground nitrogen in litter (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrnlit  (npoi)        ! global! total root litter nitrogen belowground (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totnsoi   (npoi)        ! global! total nitrogen in soil (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnmic   (npoi)        ! local! total nitrogen residing in microbial pool (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totlit    (npoi)        ! local! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totfall   (npoi)        ! local! total litterfall and root turnover (kg_C m-2/year)
    REAL(KIND=r8), INTENT(OUT  ) :: totnlit   (npoi)        ! local! total nitrogen in all litter pools (kg_N m-2)

    REAL(KIND=r8), INTENT(INOUT) :: firefac   (npoi)        ! global! factor that respresents the annual average
    REAL(KIND=r8), INTENT(INOUT) :: wtot      (npoi)        ! global! total amount of water stored in snow, soil,
    ! puddels, and on vegetation (kg_h2o)
    ! fuel dryness of a grid cell, and hence characterizes the readiness to burn

    REAL(KIND=r8), INTENT(INOUT) :: storedn (npoi)           ! global! total storage of N in soil profile (kg_N m-2) 
    REAL(KIND=r8), INTENT(INOUT) :: yrleach (npoi)           ! global! annual total amount C leached from soil profile (kg_C m-2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ynleach (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: falll   (npoi)          ! global ! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: fallr   (npoi)          ! global ! annual root litter input                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: fallw   (npoi)          ! global ! annual wood litter fall                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: clitlm  (npoi)          ! global! carbon in leaf litter pool - metabolic       (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitls  (npoi)          ! global! carbon in leaf litter pool - structural      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrm  (npoi)          ! global! carbon in fine root litter pool - metabolic  (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrs  (npoi)          ! global! carbon in fine root litter pool - structural (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwm  (npoi)          ! global! carbon in woody litter pool - metabolic      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitws  (npoi)          ! global! carbon in woody litter pool - structural     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislop(npoi)          ! global! carbon in soil - slow protected humus           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislon(npoi)          ! global! carbon in soil - slow nonprotected humus     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoipas (npoi)          ! global! carbon in soil - passive humus                   (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitll  (npoi)          ! global! carbon in leaf litter pool - lignin           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrl  (npoi)          ! global! carbon in fine root litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwl  (npoi)          ! global! carbon in woody litter pool - lignin           (kg_C m-2)


    REAL(KIND=r8), INTENT(INOUT) :: tc      (npoi)          ! global  ! coldest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: agddu   (npoi)          ! global  ! annual accumulated growing degree days for bud
    ! burst, upper canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: tempu   (npoi)          ! global  ! cold-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: agddl   (npoi)          ! global  ! annual accumulated growing degree days for bud burst,
    ! lower canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: templ   (npoi)           ! global  ! cold-phenology trigger for grasses/shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropu   (npoi)          ! global  ! drought-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropls  (npoi)          ! global  ! drought-phenology trigger for shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl4  (npoi)          ! global  ! drought-phenology trigger for c4 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl3  (npoi)          ! global  ! drought-phenology trigger for c3 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: plai    (npoi,npft)     ! global  ! total leaf area index of each plant functional type
    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iday         ! day number  (passed in)
    INTEGER, INTENT(IN   ) :: imonth         ! month number (passed in)
    INTEGER, INTENT(IN   ) :: iyear
    INTEGER, INTENT(IN   ) :: iyear0
    INTEGER, INTENT(IN   ) :: isimveg 
    INTEGER, INTENT(IN   ) :: spinmax 
    REAL(KIND=r8), INTENT(IN        ) :: ux  (npoi)
    REAL(KIND=r8), INTENT(IN        ) :: uy  (npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: taux(npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: tauy(npoi)
    REAL(KIND=r8), INTENT(INOUT) :: ts2 (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: qs2 (npoi)
    REAL(KIND=r8), INTENT(IN        ) :: deltat   (npoi)      ! absolute minimum temperature -
    ! temp on average of coldest month (C)
    REAL(KIND=r8), INTENT(INOUT) :: gdd0     (npoi)         ! growing degree days > 0C 
    REAL(KIND=r8), INTENT(INOUT) :: gdd0this (npoi)         ! annual total growing degree days for current year
    REAL(KIND=r8), INTENT(INOUT) :: tcthis   (npoi)      ! coldest monthly temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: twthis   (npoi)      ! warmest monthly temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: tcmin    (npoi)      ! coldest daily temperature of current year (C)
    REAL(KIND=r8), INTENT(INOUT) :: gdd5     (npoi)      ! growing degree days > 5C
    REAL(KIND=r8), INTENT(INOUT) :: gdd5this (npoi)      ! annual total growing degree days for current year
    REAL(KIND=r8), INTENT(IN   ) :: TminL    (npft)      ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: TminU    (npft)      ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: Twarm    (npft)      ! Temperature of warmest month (lower canopy PFTs)
    REAL(KIND=r8), INTENT(IN   ) :: GDD      (npft)      ! minimum GDD needed (base 5 C for upper canopy PFTs, 
    REAL(KIND=r8), INTENT(IN   ) :: aleaf    (npft)          ! carbon allocation fraction to leaves
    REAL(KIND=r8), INTENT(IN   ) :: awood    (npft)          ! carbon allocation fraction to wood 
    REAL(KIND=r8), INTENT(INOUT) :: cbiol    (npoi,npft) ! carbon in leaf biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(IN   ) :: aroot    (npft)         ! carbon allocation fraction to fine roots
    REAL(KIND=r8), INTENT(OUT  ) :: disturbf (npoi)         ! annual fire disturbance regime (m2/m2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: disturbo (npoi)         ! fraction of biomass pool lost every year to disturbances other than fire
    REAL(KIND=r8), INTENT(IN   ) :: specla   (npft)         ! specific leaf area (m**2/kg) 
    REAL(KIND=r8), INTENT(OUT  ) :: biomass  (npoi,npft) ! total biomass of each plant functional type  (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totlaiu  (npoi)         ! total leaf area index for the upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: totlail  (npoi)         ! total leaf area index for the lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: totbiou  (npoi)         ! total biomass in the upper canopy (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totbiol  (npoi)         ! total biomass in the lower canopy (kg_C m-2)
    REAL(KIND=r8), INTENT(IN   ) :: woodnorm                   ! value of woody biomass for upper canopy closure
    ! (ie when wood = woodnorm fu = 1.0) (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: vegtype0 (npoi)      ! annual vegetation type - ibis classification
    REAL(KIND=r8), INTENT(IN   ) :: tauwood0 (npft)      ! normal (unstressed) turnover time for wood biomass (years)
    REAL(KIND=r8), INTENT(OUT  ) :: tauwood  (npoi,npft)      ! wood biomass turnover time constant (years)
    REAL(KIND=r8), INTENT(IN   ) :: tauleaf  (npft)      ! foliar biomass turnover time constant (years)
    REAL(KIND=r8), INTENT(IN   ) :: tauroot  (npft)      ! fine root biomass turnover time constant (years)
    REAL(KIND=r8), INTENT(IN   ) :: xminlai                 ! Minimum LAI for each existing PFT
    REAL(KIND=r8), INTENT(OUT  ) :: cdisturb (npoi)         ! annual amount of vegetation carbon lost 
    ! to atmosphere due to fire  (biomass burning) (kg_C m-2/year)
    REAL(KIND=r8), INTENT(OUT  ) :: ayanpp   (npoi,npft)   ! annual above-ground npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayanpptot(npoi)             ! annual above-ground npp for ecosystem (kg-c/m**2/yr)
    !REAL(KIND=r8), INTENT(IN   ) :: garea    (npoi)   ! area of each gridcell (m**2)
    REAL(KIND=r8), INTENT(INOUT) :: tw      (npoi)      ! warmest monthly temperature (C)

    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: nstep      ! atm time step index
    !INTEGER, INTENT(IN   ) :: idayprev   ! day in month of previous timestep
    INTEGER, INTENT(IN   ) :: imonthprev ! month of previous timestep 
    !INTEGER, INTENT(IN   ) :: iyearprev  ! year of previous timestep 
    !INTEGER, INTENT(IN   ) :: idayout         ! write out daily output
    !INTEGER, INTENT(IN   ) :: imonthout  ! write out daily output
    !INTEGER, INTENT(IN   ) :: iyearout
    INTEGER, INTENT(IN   ) :: isimco2
    INTEGER, INTENT(IN   ) :: isimfire
    INTEGER, PARAMETER     :: lenMonth(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/) 
    REAL(KIND=r8), INTENT(IN   ) :: co2init
    !REAL(KIND=r8), INTENT(IN   ) :: calday             ! current julian day (1-365.99)
    !  decimals=fraction of day

    !  INTEGER :: iday         ! current day in month (1-31 or 1-30, or 1-28)
    !   integer:: imonth         ! current month (1 - 12)
    !INTEGER :: mcdate
    !REAL(KIND=r8) :: co2vmrgcm            ! modified co2 volume mixing ratio
    !REAL(KIND=r8) :: dtibis             ! time step as passed from GCM, dtime from
    INTEGER :: lenMonthly

    INTEGER :: j
    INTEGER :: i
    INTEGER :: ndyn
    terml=0.0_r8
    termu=0.0_r8
    abupd=0.0_r8
    abupi=0.0_r8
    fupdd=0.0_r8
    !
    !
    CALL Ibis    (mcsec         ,pi          ,stef         ,vonk         ,grav        , &
         tmelt        ,hfus        ,hvap         ,hsub         ,ch2o        , &
         cice         ,cair        ,cvap         ,rair         ,rvap        , &
         cappa        ,rhow        ,npoi         ,nband        ,nsoilay     , & 
         nsnolay      ,npft        ,epsilon      ,dtime        ,doalb       , &
         ginvap       ,gsuvap      ,gtrans       ,gtransu      ,gtransl     , &
         grunof       ,gdrain      ,gadjust      ,a10scalparamu,a10daylightu, &
         a10scalparaml,a10daylightl,vmax_pft     ,tau15        ,kc15        , &
         ko15         ,cimax       ,gammaub      ,alpha3       ,theta3      , &
         beta3        ,coefmub     ,coefbub      ,gsubmin      ,gammauc     , &
         coefmuc      ,coefbuc     ,gsucmin      ,gammals      ,coefmls     , & 
         coefbls      ,gslsmin     ,gammal3      ,coefml3      ,coefbl3     , &
         gsl3min      ,gammal4     ,alpha4       ,theta4       ,beta4       , &
         coefml4      ,coefbl4     ,gsl4min      ,wliqu        ,wliqumax    , & 
         wsnou        ,wsnoumax    ,tu           ,wliqs        ,wliqsmax    , & 
         wsnos        ,wsnosmax    ,ts           ,wliql        ,wliqlmax    , &  
         wsnol        ,wsnolmax    ,tl           ,topparu      ,topparl     , &
         fl           ,fu          ,lai          ,sai          ,rhoveg      , &   
         tauveg       ,orieh       ,oriev        ,wliqmin      ,wsnomin     , & 
         t12          ,tdripu      ,tblowu       ,tdrips       ,tblows      , &
         t34          ,tdripl      ,tblowl       ,ztop         ,alaiml      , &
         zbot         ,alaimu      ,froot        ,q34          ,q12         , &
         su                 ,cleaf       ,dleaf        ,ss           ,cstem       , & 
         dstem        ,sl          ,cgrass       ,ciub         ,ciuc        , &
         exist        ,csub        ,gsub         ,csuc         ,gsuc        , &
         agcub        ,agcuc       ,ancub        ,ancuc        ,totcondub   , &
         totconduc    ,cils        ,cil3         ,cil4         ,csls        , &
         gsls         ,csl3        ,gsl3         ,csl4         ,gsl4        , &
         agcls        ,agcl4       ,agcl3        ,ancls        ,ancl4       , &
         ancl3        ,totcondls   ,totcondl3    ,totcondl4    ,chu         , &
         chs          ,chl         ,frac         ,tlsub          ,z0sno       , & 
         rhos         ,consno      ,hsnotop      ,hsnomin      ,fimin       , &
         fimax        ,fi          ,tsno         ,hsno         ,sand               , &
         clay         ,poros       ,wsoi         ,wisoi        ,consoi      , &  
         zwpmax       ,wpud        ,wipud        ,wpudmax      ,qglif       , &         
         tsoi         ,hvasug      ,hvasui       ,albsav       ,albsan      , &
         tg           ,ti          ,z0soi        ,swilt        ,sfield      , &
         stressl      ,stressu     ,stresstl     ,stresstu     ,csoi        , &         
         rhosoi       ,hsoi        ,suction      ,bex          ,upsoiu      , &  
         upsoil       ,heatg       ,heati        ,hydraul      ,porosflo    , &
         ibex         ,bperm       ,hflo            ,ta           ,asurd       , &  
         asuri        ,coszen      ,solad        ,solai        ,fira               , & 
         raina        ,qa          ,psurf        ,snowa        ,ua               , &   
         o2conc       ,co2conc     ,fwetu        ,rliqu        ,fwets       , &
         rliqs        ,fwetl       ,rliql        ,solu        , &         
         sols         ,soll        ,solg         ,soli         ,scalcoefl   , &
         scalcoefu    ,indsol      ,albsod       ,albsoi       ,albsnd      , &
         albsni       ,relod       ,reloi        ,reupd        ,reupi       , &
         ablod        ,abloi              ,flodd        ,dummy        ,flodi       , &
         floii         ,terml       ,termu        ,abupd        ,abupi       , & 
         fupdd        ,fupdi       ,fupii        ,sol2d        ,sol2i       , & 
         sol3d        ,sol3i       ,firb         ,firs         ,firu               , &
         firl         ,firg        ,firi         ,snowg        ,tsnowg      , &    
         tsnowl       ,pfluxl      ,raing        ,traing       ,trainl      , &
         snowl        ,tsnowu      ,pfluxu       ,rainu        ,trainu      , &   
         snowu        ,pfluxs      ,rainl        ,tfac         ,cp          , &  
         za                 ,bdl         ,dil          ,z3           ,z4          , & 
         z34                 ,exphl       ,expl         ,displ        ,bdu         , &
         diu          ,z1          ,z2           ,z12          ,exphu       , &
         expu         ,dispu       ,alogg        ,alogi        ,alogav      , & 
         alog4        ,alog3       ,alog2        ,alog1        ,aloga       , &          
         u2           ,alogu       ,alogl        ,richl        ,straml      , &
         strahl       ,richu       ,stramu       ,strahu       ,u1          , & 
         u12          ,u3          ,u34          ,u4           ,cu          , &
         cl                 ,sg          ,si           ,fwetux       ,fwetsx      , &  
         fwetlx       ,fsena       ,fseng        ,fseni        ,fsenu       , &
         fsens        ,fsenl       ,fvapa        ,fvaput       ,fvaps       , &
         fvaplw       ,fvaplt      ,fvapg        ,fvapi        ,fvapuw      , &
         td                 ,vzero       ,ndaypy       ,nppdummy     ,tgpp        , & 
         tgpptot      ,tnpp        ,cbiow        ,sapfrac      ,cbior       , & 
         tnpptot      ,tco2root    ,tneetot      ,tco2mic      ,a10td       , &
         a10ancub     ,a10ancuc    ,a10ancls     ,a10ancl3     ,a10ancl4    , & 
         ndtimes      ,adrain       ,adsnow      , &
         adaet        ,adtrunoff   ,adsrunoff    ,addrainage   ,adrh        , &
         adsnod       ,adsnof      ,adwsoi       ,adtsoi       ,adwisoi     , &
         adtlaysoi    ,adwlaysoi   ,adwsoic      ,adtsoic      ,adco2mic    , &
         adco2root    ,adco2soi    ,adco2ratio   ,adnmintot    ,decompl     , &    
         decomps      ,tnmin       ,ndaypm       ,nmtimes     , & 
         amrain       ,amsnow      ,amaet        ,amtrunoff    ,amsrunoff   , &
         amdrainage   ,amtemp      ,amqa         , &
         amsolar      ,amirup      ,amirdown     ,amsens       ,amlatent    , &  
         amlaiu       ,amlail      ,amtsoi       ,amwsoi       ,amwisoi     , &   
         amvwc        ,amawc       ,amsnod       ,amsnof       ,amnpp       , &
         amnpptot     ,amco2mic    ,amco2root    ,amco2soi     ,amco2ratio  , &
         amneetot     ,amnmintot   ,nytimes      ,aysolar      ,ayirup      , &
         ayirdown     ,aysens      ,aylatent     ,ayprcp       ,ayaet       , &  
         aytrans      ,aytrunoff   ,aysrunoff    ,aydrainage   ,aydwtot     , & 
         aywsoi       ,aywisoi     ,aytsoi       ,ayvwc        ,ayawc       , &  
         aystresstu   ,aystresstl  ,aygpp            ,aygpptot     ,aynpp       , & 
         aynpptot     ,ayco2mic    ,ayco2root    ,ayco2soi     ,ayneetot    , &
         ayrootbio    ,aynmintot   ,ayalit       ,ayblit       ,aycsoi      , & 
         aycmic       ,ayanlit     ,aybnlit      ,aynsoi       ,ayalbedo    , &
         totalit     , &
         totrlit      ,totcsoi     ,totcmic      ,totanlit     ,totrnlit    , &
         totnsoi      ,totnmic     ,totlit       ,totfall      ,totnlit     , &
         firefac      ,wtot        ,storedn      ,yrleach      ,ynleach     , & 
         falll        ,fallr       ,fallw        ,clitlm       ,clitls      , &
         clitrm       ,clitrs      ,clitwm       ,clitws       ,csoislop    , &
         csoislon     ,csoipas     ,clitll       ,clitrl       ,clitwl      , &  
         tc           ,agddu              ,tempu            ,agddl        ,templ       , &
         dropu        ,dropls      ,dropl4       ,dropl3       ,plai               , &
         iday        ,imonth       ,iyear        ,iyear0      , &
         isimveg      ,spinmax     ,amts2        ,amtransu    ,amtransl     , &
         amsuvap      ,aminvap     ,amalbedo     ,amtsoil     ,amwsoil      , & 
         amwisoil     ,ux          ,uy           ,taux        ,tauy         , &
         ts2                 ,qs2         ,gdd0this     ,gdd5this    ,bstar  )


    DO i = 1, npoi
       z0    (i) =  EXP(((alogu(i))))
       ustar (i) =  ua (i) * cu (i)
       hc    (i) =  (fsena  (i) + fsenu    (i) + fsens     (i) + fsenl    (i) +fsenl(i) )*dtime
       hg    (i) =  (fseng  (i) + fseni    (i)                                          )*dtime

       !hltm   = 2.52e6_r8            !  latent heat of vaporization (J kg^-1)

       ec    (i) =  (fvaput (i) + fvaps    (i) + fvaplw    (i)  + fvaplt (i) +fvapuw (i))*hltm*dtime
       eg    (i) =  (fvapg  (i) + fvapi    (i)                                                 )*hltm*dtime
    END DO


    IF (nstep > 20) THEN
       ! IF (nstep /= 0) THEN
       !        
       lenMonthly=lenMonth(imonth)
       IF(MOD(REAL(iyear,KIND=r8),4.0_r8) == 0.0_r8 .AND. imonth == 2) lenMonthly=29
       IF ((mcsec == 86400.0_r8-(dtime/2.0_r8)) .AND. (iday == lenMonthly) ) THEN 
          !
          !     end of calculations done once a month
          !
          DO i = 1, npoi
             tcthis(i) = MIN (tcthis(i), (amts2(i) - 273.160_r8))
             twthis(i) = MAX (twthis(i), (amts2(i) - 273.160_r8))
          END DO
       END IF
       !
       ! calculations done once a year
       !
       IF ((nstep > 360) .AND.(mcsec == 86400.0_r8-(dtime/2.0_r8)) .AND. (iday == lenMonthly) .AND. (imonth == 12) ) THEN
          !         IF ((mcsec == 0.0_r8) .and. (iday == 1) .and. (imonth == 1) ) THEN

          !
          ! get new o2 and co2 concentrations for this year
          !
          IF (isimco2.EQ.1) CALL co2 (co2init,  &! INTENT(IN   )
               co2conc,  &! INTENT(OUT  )
               iyear     )! INTENT(IN   )

          !
          ! perform vegetation dynamics
          !
          ndyn = 1
          DO j = 1, ndyn

             IF (isimveg /= 0) CALL dynaveg2 (isimfire , &! INTENT(IN   )        dynaveg1 (isimfire , &! INTENT(IN   )
                  tauwood0  , &! INTENT(IN   )                  tauwood0 , &! INTENT(IN   )
                  tauwood   , &! INTENT(OUT  )                  tauwood  , &! INTENT(OUT  )
                  tauleaf   , &! INTENT(IN   )                  tauleaf  , &! INTENT(IN   )
                  tauroot   , &! INTENT(IN   )                  tauroot  , &! INTENT(IN   )
                  xminlai   , &! INTENT(IN   )                  xminlai  , &! INTENT(IN   )
                  falll     , &! INTENT(OUT  )                  falll    , &! INTENT(OUT  )
                  fallr     , &! INTENT(OUT  )                  fallr    , &! INTENT(OUT  )
                  fallw     , &! INTENT(OUT  )                  fallw    , &! INTENT(OUT  )
                  cdisturb  , &! INTENT(OUT  )                  cdisturb , &! INTENT(OUT  )
                  exist     , &! INTENT(OUT  )                  exist    , &! INTENT(IN   )
                  aleaf     , &! INTENT(IN   )                  aleaf    , &! INTENT(IN   )
                  awood     , &! INTENT(IN   )                  awood    , &! INTENT(IN   )
                  cbiol     , &! INTENT(INOUT) global           cbiol    , &! INTENT(INOUT) global
                  cbior     , &! INTENT(INOUT) global           cbior    , &! INTENT(INOUT) global
                  cbiow     , &! INTENT(INOUT) global           cbiow    , &! INTENT(INOUT) global
                  aroot     , &! INTENT(IN   )                  aroot    , &! INTENT(IN   )
                  disturbf  , &! INTENT(OUT  )                  disturbf , &! INTENT(OUT  )
                  disturbo  , &! INTENT(OUT  )                  disturbo , &! INTENT(OUT  )
                  firefac   , &! INTENT(IN   )                  firefac  , &! INTENT(IN   )
                  totlit    , &! INTENT(IN   )                  totlit   , &! INTENT(IN   )
                  specla    , &! INTENT(IN   )                  specla   , &! INTENT(IN   )
                  plai      , &! INTENT(INOUT) local                 plai          , &! INTENT(INOUT) local
                  biomass   , &! INTENT(OUT  )                  biomass  , &! INTENT(OUT  )
                  totlaiu   , &! INTENT(INOUT) local                 totlaiu  , &! INTENT(INOUT) local
                  totlail   , &! INTENT(INOUT) local                 totlail  , &! INTENT(INOUT) local
                  totbiou   , &! INTENT(INOUT) local                 totbiou  , &! INTENT(INOUT) local
                  totbiol   , &! INTENT(OUT  )                  totbiol  , &! INTENT(OUT  )
                  fu             , &! INTENT(OUT  )                  fu          , &! INTENT(OUT  )
                  woodnorm  , &! INTENT(IN   )                  woodnorm , &! INTENT(IN   )
                  fl             , &! INTENT(OUT  )                  fl          , &! INTENT(OUT  )
                  zbot      , &! INTENT(OUT  )                  zbot          , &! INTENT(OUT  )
                  ztop      , &! INTENT(OUT  )                  ztop          , &! INTENT(OUT  )
                  sai       , &! INTENT(OUT  )                  sai          , &! INTENT(OUT  )
                  sapfrac   , &! INTENT(OUT  )                  sapfrac  , &! INTENT(OUT  )
                  vegtype0  , &! INTENT(OUT  )                  vegtype0 , &! INTENT(OUT  )
                  gdd5      , &! INTENT(IN   )                  gdd5          , &! INTENT(IN   )
                  gdd0      , &! INTENT(IN   )                  gdd0          , &! INTENT(IN   )
                  aynpp     , &! INTENT(INOUT) global           aynpp    , &! INTENT(INOUT) global
                  ayanpp    , &! INTENT(OUT  )                  ayanpp   , &! INTENT(OUT  )
                  ayneetot  , &! INTENT(INOUT) global           ayneetot , &! INTENT(INOUT) global
                  ayanpptot , &! INTENT(OUT  )                  ayanpptot, &! INTENT(OUT  )
                  aynpptot  , &! INTENT(OUT  )                  npoi          , &!
                  ayco2mic  , &! INTENT(IN   )                  npft            )! , isim_ac, year)
                  npoi      , &!
                  npft        )! , isim_ac, year)!

             !IF (isimveg /= 0) CALL dynaveg1 (isimfire , &! INTENT(IN        )
             !                             tauwood0 , &! INTENT(IN   )
             !                             tauwood  , &! INTENT(OUT  )
             !                             tauleaf  , &! INTENT(IN   )
             !                             tauroot  , &! INTENT(IN   )
             !                             xminlai  , &! INTENT(IN   )
             !                             falll    , &! INTENT(OUT  )
             !                             fallr    , &! INTENT(OUT  )
             !                             fallw    , &! INTENT(OUT  )
             !                             cdisturb , &! INTENT(OUT  )
             !                             exist    , &! INTENT(IN   )
             !                             aleaf    , &! INTENT(IN   )
             !                             awood    , &! INTENT(IN   )
             !                             cbiol    , &! INTENT(INOUT) global
             !                             cbior    , &! INTENT(INOUT) global
             !                             cbiow    , &! INTENT(INOUT) global
             !                             aroot    , &! INTENT(IN   )
             !                             disturbf , &! INTENT(OUT  )
             !                             disturbo , &! INTENT(OUT  )
             !                             firefac  , &! INTENT(IN   )
             !                             totlit   , &! INTENT(IN   )
             !                             specla   , &! INTENT(IN   )
             !                             plai     , &! INTENT(INOUT) local
             !                             biomass  , &! INTENT(OUT  )
             !                             totlaiu  , &! INTENT(INOUT) local
             !                             totlail  , &! INTENT(INOUT) local
             !                             totbiou  , &! INTENT(INOUT) local
             !                             totbiol  , &! INTENT(OUT  )
             !                             fu       , &! INTENT(OUT  )
             !                             woodnorm , &! INTENT(IN   )
             !                             fl       , &! INTENT(OUT  )
             !                             zbot     , &! INTENT(OUT  )
             !                             ztop     , &! INTENT(OUT  )
             !                             sai      , &! INTENT(OUT  )
             !                             sapfrac  , &! INTENT(OUT  )
             !                             vegtype0 , &! INTENT(OUT  )
             !                             gdd5     , &! INTENT(IN   )
             !                             gdd0     , &! INTENT(IN   )
             !                             aynpp    , &! INTENT(INOUT) global
             !                             ayanpp   , &! INTENT(OUT  )
             !                             ayneetot , &! INTENT(INOUT) global
             !                             ayanpptot, &! INTENT(OUT  )
             !                             npoi     , &!
             !                             npft        )! , isim_ac, year)
             !
             !
          END DO

          !
          !
          !     recalculate bioclimatic parameters (used in dynaveg, calculated
          !     even in fixed vegetation case when fixed vegetation is an
          !     initialisation of a dynamic run)
          !
          CALL climanl2(TminL    , &! INTENT(IN   )
               TminU    , &! INTENT(IN   )
               Twarm    , &! INTENT(IN   )
               GDD      , &! INTENT(IN   )
               gdd0     , &! INTENT(INOUT)
               gdd0this , &! INTENT(IN   )
               tc       , &! INTENT(INOUT)
               tw       , &! INTENT(INOUT)
               tcthis   , &! INTENT(IN   )
               twthis   , &! INTENT(IN   )
               tcmin    , &! INTENT(INOUT) local
               gdd5     , &! INTENT(INOUT) local
               gdd5this , &! INTENT(IN   )
               exist    , &! INTENT(INOUT)
               deltat   , &! INTENT(IN   )
               npoi     , &! INTENT(IN   )
               npft       )! INTENT(IN   )

       END IF

       IF (imonthprev .NE. imonth) THEN
          !
          ! write restart files
          !
          !CALL wrestart (mdcur, imonthprev, iyearprev, iyear0)
       END IF
       !
       !     End of test on 1st time step
       !
    END IF
  END SUBROUTINE IbisDrv






  SUBROUTINE Ibis   (mcsec        ,pi          ,stef         ,vonk         ,grav        , &
       tmelt        ,hfus        ,hvap         ,hsub         ,ch2o        , &
       cice         ,cair        ,cvap         ,rair         ,rvap        , &
       cappa        ,rhow        ,npoi         ,nband        ,nsoilay     , & 
       nsnolay      ,npft        ,epsilon      ,dtime        ,doalb       , &
       ginvap       ,gsuvap      ,gtrans       ,gtransu      ,gtransl     , &
       grunof       ,gdrain      ,gadjust      ,a10scalparamu,a10daylightu, &
       a10scalparaml,a10daylightl,vmax_pft     ,tau15        ,kc15        , &
       ko15         ,cimax       ,gammaub      ,alpha3       ,theta3      , &
       beta3        ,coefmub     ,coefbub      ,gsubmin      ,gammauc     , &
       coefmuc      ,coefbuc     ,gsucmin      ,gammals      ,coefmls     , & 
       coefbls      ,gslsmin     ,gammal3      ,coefml3      ,coefbl3     , &
       gsl3min      ,gammal4     ,alpha4       ,theta4       ,beta4       , &
       coefml4      ,coefbl4     ,gsl4min      ,wliqu        ,wliqumax    , & 
       wsnou        ,wsnoumax    ,tu           ,wliqs        ,wliqsmax    , & 
       wsnos        ,wsnosmax    ,ts           ,wliql        ,wliqlmax    , &  
       wsnol        ,wsnolmax    ,tl           ,topparu      ,topparl     , &
       fl           ,fu          ,lai          ,sai          ,rhoveg      , &   
       tauveg       ,orieh       ,oriev        ,wliqmin      ,wsnomin     , & 
       t12          ,tdripu      ,tblowu       ,tdrips       ,tblows      , &
       t34          ,tdripl      ,tblowl       ,ztop         ,alaiml      , &
       zbot         ,alaimu      ,froot        ,q34          ,q12         , &
       su                 ,cleaf       ,dleaf        ,ss           ,cstem       , & 
       dstem        ,sl          ,cgrass       ,ciub         ,ciuc        , &
       exist        ,csub        ,gsub         ,csuc         ,gsuc        , &
       agcub        ,agcuc       ,ancub        ,ancuc        ,totcondub   , &
       totconduc    ,cils        ,cil3         ,cil4         ,csls        , &
       gsls         ,csl3        ,gsl3         ,csl4         ,gsl4        , &
       agcls        ,agcl4       ,agcl3        ,ancls        ,ancl4       , &
       ancl3       ,totcondls    ,totcondl3    ,totcondl4    ,chu         , &
       chs          ,chl         ,frac         ,tlsub          ,z0sno       , & 
       rhos         ,consno      ,hsnotop      ,hsnomin      ,fimin       , &
       fimax        ,fi          ,tsno         ,hsno         ,sand               , &
       clay         ,poros       ,wsoi         ,wisoi        ,consoi      , &  
       zwpmax       ,wpud        ,wipud        ,wpudmax      ,qglif       , &         
       tsoi         ,hvasug      ,hvasui       ,albsav       ,albsan      , &
       tg           ,ti          ,z0soi        ,swilt        ,sfield      , &
       stressl      ,stressu     ,stresstl     ,stresstu     ,csoi        , &         
       rhosoi       ,hsoi        ,suction      ,bex          ,upsoiu      , &  
       upsoil       ,heatg       ,heati        ,hydraul      ,porosflo    , &
       ibex         ,bperm       ,hflo            ,ta           ,asurd       , &  
       asuri        ,coszen      ,solad        ,solai        ,fira               , & 
       raina        ,qa          ,psurf        ,snowa        ,ua               , &   
       o2conc       ,co2conc     ,fwetu        ,rliqu        ,fwets       , &
       rliqs        ,fwetl       ,rliql        ,solu        , &         
       sols         ,soll        ,solg         ,soli         ,scalcoefl   , &
       scalcoefu    ,indsol      ,albsod       ,albsoi       ,albsnd      , &
       albsni       ,relod       ,reloi        ,reupd        ,reupi       , &
       ablod        ,abloi              ,flodd        ,dummy        ,flodi       , &
       floii         ,terml       ,termu        ,abupd        ,abupi       , & 
       fupdd        ,fupdi       ,fupii        ,sol2d        ,sol2i       , & 
       sol3d        ,sol3i       ,firb         ,firs         ,firu               , &
       firl         ,firg        ,firi         ,snowg        ,tsnowg      , &    
       tsnowl       ,pfluxl      ,raing        ,traing       ,trainl      , &
       snowl        ,tsnowu      ,pfluxu       ,rainu        ,trainu      , &   
       snowu        ,pfluxs      ,rainl        ,tfac         ,cp          , &  
       za                 ,bdl         ,dil          ,z3           ,z4          , & 
       z34                 ,exphl       ,expl         ,displ        ,bdu         , &
       diu          ,z1          ,z2           ,z12          ,exphu       , &
       expu         ,dispu       ,alogg        ,alogi        ,alogav      , & 
       alog4        ,alog3       ,alog2        ,alog1        ,aloga       , &          
       u2           ,alogu       ,alogl        ,richl        ,straml      , &
       strahl       ,richu       ,stramu       ,strahu       ,u1          , & 
       u12          ,u3          ,u34          ,u4           ,cu          , &
       cl                 ,sg          ,si           ,fwetux       ,fwetsx      , &  
       fwetlx       ,fsena       ,fseng        ,fseni        ,fsenu       , &
       fsens        ,fsenl       ,fvapa        ,fvaput       ,fvaps       , &
       fvaplw       ,fvaplt      ,fvapg        ,fvapi        ,fvapuw      , &
       td                 ,vzero       ,ndaypy       ,nppdummy     ,tgpp        , & 
       tgpptot      ,tnpp        ,cbiow        ,sapfrac      ,cbior       , & 
       tnpptot      ,tco2root    ,tneetot      ,tco2mic      ,a10td       , &
       a10ancub     ,a10ancuc    ,a10ancls     ,a10ancl3     ,a10ancl4    , & 
       ndtimes      ,adrain       ,adsnow      , &
       adaet        ,adtrunoff   ,adsrunoff    ,addrainage   ,adrh        , &
       adsnod       ,adsnof      ,adwsoi       ,adtsoi       ,adwisoi     , &
       adtlaysoi    ,adwlaysoi   ,adwsoic      ,adtsoic      ,adco2mic    , &
       adco2root    ,adco2soi    ,adco2ratio   ,adnmintot    ,decompl     , &    
       decomps      ,tnmin       ,ndaypm       ,nmtimes     , & 
       amrain       ,amsnow      ,amaet        ,amtrunoff    ,amsrunoff   , &
       amdrainage   ,amtemp      ,amqa         ,&
       amsolar      ,amirup      ,amirdown     ,amsens       ,amlatent    , &  
       amlaiu       ,amlail      ,amtsoi       ,amwsoi       ,amwisoi     , &   
       amvwc        ,amawc       ,amsnod       ,amsnof       ,amnpp       , &
       amnpptot     ,amco2mic    ,amco2root    ,amco2soi     ,amco2ratio  , &
       amneetot     ,amnmintot   ,nytimes      ,aysolar      ,ayirup      , &
       ayirdown     ,aysens      ,aylatent     ,ayprcp       ,ayaet       , &  
       aytrans      ,aytrunoff   ,aysrunoff    ,aydrainage   ,aydwtot     , & 
       aywsoi       ,aywisoi     ,aytsoi       ,ayvwc        ,ayawc       , &  
       aystresstu   ,aystresstl  ,aygpp            ,aygpptot     ,aynpp       , & 
       aynpptot     ,ayco2mic    ,ayco2root    ,ayco2soi     ,ayneetot    , &
       ayrootbio    ,aynmintot   ,ayalit       ,ayblit       ,aycsoi      , & 
       aycmic       ,ayanlit     ,aybnlit      ,aynsoi       ,ayalbedo    , &
       totalit      , &
       totrlit      ,totcsoi     ,totcmic      ,totanlit     ,totrnlit    , &
       totnsoi      ,totnmic     ,totlit       ,totfall      ,totnlit     , &
       firefac      ,wtot        ,storedn      ,yrleach      ,ynleach     , & 
       falll        ,fallr       ,fallw        ,clitlm       ,clitls      , &
       clitrm       ,clitrs      ,clitwm       ,clitws       ,csoislop    , &
       csoislon     ,csoipas     ,clitll       ,clitrl       ,clitwl      , &  
       tc           ,agddu              ,tempu            ,agddl        ,templ       , &
       dropu        ,dropls      ,dropl4       ,dropl3       ,plai               , &
       iday        ,imonth       ,iyear        ,iyear0       , &
       isimveg      ,spinmax     ,amts2        ,amtransu     ,amtransl    , &
       amsuvap      ,aminvap     ,amalbedo     ,amtsoil      ,amwsoil     , & 
       amwisoil     , ux              ,uy           ,taux         ,tauy        , &
       ts2                 ,qs2         ,gdd0this     ,gdd5this     ,bstar)

    IMPLICIT NONE
    REAL(KIND=r8), INTENT(IN        ) :: mcsec!global  ! current seconds in day (0 - (86400 - dtime))
    REAL(KIND=r8), INTENT(IN        ) :: pi   !global
    REAL(KIND=r8), INTENT(IN        ) :: stef !global  ! stefan-boltzmann constant (W m-2 K-4)
    REAL(KIND=r8), INTENT(IN        ) :: vonk !global  ! von karman constant (dimensionless)
    REAL(KIND=r8), INTENT(IN        ) :: grav !global  ! gravitational acceleration (m s-2)
    REAL(KIND=r8), INTENT(IN        ) :: tmelt!global  ! freezing point of water (K)
    REAL(KIND=r8), INTENT(IN        ) :: hfus !global  ! latent heat of fusion of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hvap !global  ! latent heat of vaporization of water (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: hsub !global  ! latent heat of sublimation of ice (J kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: ch2o !global  ! specific heat of liquid water (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cice !global  ! specific heat of ice (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cair !global  ! specific heat of dry air at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cvap !global  ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rair !global  ! gas constant for dry air (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: rvap !global  ! gas constant for water vapor (J deg-1 kg-1)
    REAL(KIND=r8), INTENT(IN        ) :: cappa!global  ! rair/cair
    REAL(KIND=r8), INTENT(IN        ) :: rhow !global  ! density of liquid water (all types) (kg m-3)

    ! 
    !
    INTEGER, INTENT(IN   ) :: npoi   !global  
    INTEGER, INTENT(IN   ) :: nband  !global  
    INTEGER, INTENT(IN   ) :: nsoilay!global   ! number of soil layers
    INTEGER, INTENT(IN   ) :: nsnolay!global   ! number of snow layers
    INTEGER, INTENT(IN   ) :: npft   !global   ! number of plant functional types
    REAL(KIND=r8), INTENT(IN   ) :: epsilon!global   ! small quantity to avoid zero-divides and other
    ! truncation or machine-limit troubles with small
    ! values. should be slightly greater than o(1)
    ! machine precision
    REAL(KIND=r8), INTENT(IN   )  :: dtime !global   ! model timestep (seconds)
    LOGICAL, INTENT(IN   )  :: doalb !global    ! true if surface albedo calculation time step

    !      INCLUDE 'comhyd.h'
    REAL(KIND=r8), INTENT(OUT  ) :: ginvap (npoi)!local ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gsuvap (npoi)!local ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gtrans (npoi)!local ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransu(npoi)!local ! transpiration from upper canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(OUT  ) :: gtransl(npoi)!local ! transpiration from lower canopy (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: grunof (npoi)!local ! surface runoff rate (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gdrain (npoi)!local ! drainage rate out of bottom of lowest soil layer (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: gadjust(npoi)!local ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
    !      INCLUDE 'comsum.h'
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparamu(npoi)!global ! 10-day average day-time scaling parameter - upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightu (npoi)!global ! 10-day average day-time PAR - upper canopy (micro-Ein m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10scalparaml(npoi)!global ! 10-day average day-time scaling parameter - lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: a10daylightl (npoi)!global ! 10-day average day-time PAR - lower canopy (micro-Ein m-2 s-1)
    !      INCLUDE 'compft.h'
    REAL(KIND=r8), INTENT(IN   ) :: vmax_pft(npft)!global ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]
    REAL(KIND=r8), INTENT(IN   ) :: tau15           !global ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: kc15           !global ! co2 kinetic parameter (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: ko15           !global ! o2 kinetic parameter (mol/mol) 
    REAL(KIND=r8), INTENT(IN   ) :: cimax           !global ! maximum value for ci (needed for model stability)
    REAL(KIND=r8), INTENT(IN   ) :: gammaub           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha3           !global ! intrinsic quantum efficiency for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: beta3           !global ! photosynthesis coupling coefficient for C3 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefmub           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbub           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsubmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammauc           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmuc           !global ! 'm' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: coefbuc           !global ! 'b' coefficient for stomatal conductance relationship  
    REAL(KIND=r8), INTENT(IN   ) :: gsucmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammals           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefmls           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbls           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gslsmin           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal3           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: coefml3           !global ! 'm' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: coefbl3           !global ! 'b' coefficient for stomatal conductance relationship 
    REAL(KIND=r8), INTENT(IN   ) :: gsl3min           !global ! absolute minimum stomatal conductance
    REAL(KIND=r8), INTENT(IN   ) :: gammal4           !global ! leaf respiration coefficient
    REAL(KIND=r8), INTENT(IN   ) :: alpha4           !global ! intrinsic quantum efficiency for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: theta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless) 
    REAL(KIND=r8), INTENT(IN   ) :: beta4           !global ! photosynthesis coupling coefficient for C4 plants (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: coefml4           !global ! 'm' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: coefbl4           !global ! 'b' coefficient for stomatal conductance relationship
    REAL(KIND=r8), INTENT(IN   ) :: gsl4min           !global ! absolute minimum stomatal conductance
    !      include 'comveg.h'
    REAL(KIND=r8), INTENT(INOUT) :: wliqu    (npoi)!global ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqumax              !global ! maximum intercepted water on a unit upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnou    (npoi)!global ! intercepted frozen h2o (snow) on upper canopy leaf area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnoumax              !global ! intercepted snow capacity for upper canopy leaves (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tu       (npoi)!global ! temperature of upper canopy leaves (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliqs    (npoi)!global ! intercepted liquid h2o on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqsmax              !global ! maximum intercepted water on a unit upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnos    (npoi)!global ! intercepted frozen h2o (snow) on upper canopy stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnosmax              !global ! intercepted snow capacity for upper canopy stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: ts       (npoi)!global ! temperature of upper canopy stems (K)
    REAL(KIND=r8), INTENT(INOUT) :: wliql    (npoi)!global ! intercepted liquid h2o on lower canopy leaf and stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wliqlmax              !global ! maximum intercepted water on a unit lower canopy stem & leaf area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnol    (npoi)!global ! intercepted frozen h2o (snow) on lower canopy leaf & stem area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wsnolmax              !global ! intercepted snow capacity for lower canopy leaves & stems (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: tl       (npoi)!global ! temperature of lower canopy leaves & stems(K)
    REAL(KIND=r8), INTENT(INOUT) :: topparu  (npoi)!local  ! total photosynthetically active raditaion absorbed 
    ! by top leaves of upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: topparl  (npoi)!local  ! total photosynthetically active raditaion absorbed
    ! by top leaves of lower canopy (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fl       (npoi)!global ! fraction of snow-free area covered by lower  canopy
    REAL(KIND=r8), INTENT(IN   ) :: fu       (npoi)!global ! fraction of overall area covered by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: lai      (npoi,2)!global ! canopy single-sided leaf area index (area leaf/area veg)
    REAL(KIND=r8), INTENT(IN   ) :: sai      (npoi,2)!global ! current single-sided stem area index
    REAL(KIND=r8), INTENT(IN   ) :: rhoveg   (nband,2)!global ! reflectance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: tauveg   (nband,2)  ! transmittance of an average leaf/stem
    REAL(KIND=r8), INTENT(IN   ) :: orieh    (2) ! fraction of leaf/stems with horizontal orientation
    REAL(KIND=r8), INTENT(IN   ) :: oriev    (2) ! fraction of leaf/stems with vertical
    REAL(KIND=r8), INTENT(INOUT) :: wliqmin      ! local ! minimum intercepted water on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wsnomin      ! local ! minimum intercepted snow on unit vegetated area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: t12      (npoi) !global ! air temperature at z12 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripu       ! global ! decay time for dripoff of liquid intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowu          ! global ! decay time for blowoff of snow intercepted by upper canopy leaves (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tdrips          ! global ! decay time for dripoff of liquid intercepted by upper canopy stems (sec) 
    REAL(KIND=r8), INTENT(IN   ) :: tblows          ! global ! decay time for blowoff of snow intercepted by upper canopy stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: t34      (npoi)! global ! air temperature at z34 (K)
    REAL(KIND=r8), INTENT(IN   ) :: tdripl       ! global ! decay time for dripoff of liquid intercepted
    ! by lower canopy leaves & stem (sec)
    REAL(KIND=r8), INTENT(IN   ) :: tblowl          ! global          ! decay time for blowoff of snow intercepted by lower canopy leaves & stems (sec)
    REAL(KIND=r8), INTENT(INOUT) :: ztop     (npoi,2) ! global  ! height of plant top above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaiml          ! global ! lower canopy leaf & stem maximum area (2 sided) for
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: zbot     (npoi,2) ! global  ! height of lowest branches above ground (m)
    REAL(KIND=r8), INTENT(IN   ) :: alaimu                ! global  ! upper canopy leaf & stem area (2 sided) for 
    ! normalization of drag coefficient (m2 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: froot    (nsoilay,2)! global! fraction of root in soil layer 
    REAL(KIND=r8), INTENT(INOUT) :: q34      (npoi)         ! global! specific humidity of air at z34
    REAL(KIND=r8), INTENT(INOUT) :: q12      (npoi)         ! global! specific humidity of air at z12
    REAL(KIND=r8), INTENT(INOUT) :: su       (npoi)         ! local ! air-vegetation transfer coefficients (*rhoa) for
    ! upper canopy leaves (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cleaf                 ! global ! empirical constant in upper canopy leaf-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dleaf    (2)         ! global ! typical linear leaf dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: ss       (npoi)         ! local! air-vegetation transfer coefficients (*rhoa) for 
    ! upper canopy stems (m s-1 * kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cstem                 ! global ! empirical constant in upper canopy stem-air 
    ! aerodynamic transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(IN   ) :: dstem    (2)         ! global ! typical linear stem dimension in aerodynamic transfer coefficient (m)
    REAL(KIND=r8), INTENT(INOUT) :: sl       (npoi)         ! local! air-vegetation transfer coefficients (*rhoa) for 
    ! lower canopy leaves & stems (m s-1*kg m-3) (A39a Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(IN   ) :: cgrass                 ! global ! empirical constant in lower canopy-air aerodynamic 
    ! transfer coefficient (m s-0.5) (A39a Pollard & Thompson 95)
    REAL(KIND=r8), INTENT(INOUT) :: ciub     (npoi)         ! global ! intercellular co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: ciuc     (npoi)         ! global ! intercellular co2 concentration - conifer        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(IN   ) :: exist    (npoi,npft)! global ! probability of existence of each plant functional type in a gridcell
    REAL(KIND=r8), INTENT(INOUT) :: csub     (npoi)         ! global ! leaf boundary layer co2 concentration - broadleaf (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsub     (npoi)         ! global ! upper canopy stomatal conductance - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csuc     (npoi)         ! global ! leaf boundary layer co2 concentration - conifer   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsuc     (npoi)         ! global ! upper canopy stomatal conductance - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcub    (npoi)         ! local  ! canopy average gross photosynthesis rate - broadleaf  (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcuc    (npoi)         ! local  ! canopy average gross photosynthesis rate - conifer    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancub    (npoi)         ! local ! canopy average net photosynthesis rate - broadleaf    (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancuc    (npoi)         ! local ! canopy average net photosynthesis rate - conifer          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondub(npoi)         ! local ! 
    REAL(KIND=r8), INTENT(INOUT) :: totconduc(npoi)         ! local !
    REAL(KIND=r8), INTENT(INOUT) :: cils     (npoi)         ! global ! intercellular co2 concentration - shrubs        (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil3     (npoi)         ! global ! intercellular co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: cil4     (npoi)         ! global ! intercellular co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: csls     (npoi)         ! global ! leaf boundary layer co2 concentration - shrubs   (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsls     (npoi)         ! global ! lower canopy stomatal conductance - shrubs     (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl3     (npoi)         ! global ! leaf boundary layer co2 concentration - c3 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl3     (npoi)         ! global ! lower canopy stomatal conductance - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: csl4     (npoi)         ! global ! leaf boundary layer co2 concentration - c4 plants (mol_co2/mol_air)
    REAL(KIND=r8), INTENT(INOUT) :: gsl4     (npoi)         ! global ! lower canopy stomatal conductance - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcls    (npoi)         ! local  ! canopy average gross photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl4    (npoi)         ! local ! canopy average gross photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: agcl3    (npoi)         ! local ! canopy average gross photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancls    (npoi)         ! local ! canopy average net photosynthesis rate - shrubs          (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl4    (npoi)         ! local ! canopy average net photosynthesis rate - c4 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: ancl3    (npoi)         ! local ! canopy average net photosynthesis rate - c3 grasses   (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: totcondls(npoi)         ! local ! 
    REAL(KIND=r8), INTENT(INOUT) :: totcondl3(npoi)         ! local !
    REAL(KIND=r8), INTENT(INOUT) :: totcondl4(npoi)         ! local !
    REAL(KIND=r8), INTENT(IN   ) :: chu                 ! global ! heat capacity of upper canopy leaves per unit leaf area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chs                 ! global ! heat capacity of upper canopy stems per unit stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(IN   ) :: chl                 ! global ! heat capacity of lower canopy leaves & stems per unit leaf/stem area (J kg-1 m-2)
    REAL(KIND=r8), INTENT(INOUT) :: frac     (npoi,npft)! global ! fraction of canopy occupied by each plant functional type
    REAL(KIND=r8), INTENT(INOUT) :: tlsub    (npoi)         ! global ! temperature of lower canopy vegetation buried by snow (K)
    !      INCLUDE 'comsat.h'    
    !      include 'comsno.h'
    REAL(KIND=r8), INTENT(IN   ) :: z0sno  ! global ! roughness length of snow surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: rhos   ! global ! density of snow (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: consno ! global ! thermal conductivity of snow (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: hsnotop! global ! thickness of top snow layer (m)
    REAL(KIND=r8), INTENT(IN   ) :: hsnomin! global ! minimum total thickness of snow (m)
    REAL(KIND=r8), INTENT(IN   ) :: fimin  ! global ! minimum fractional snow cover
    REAL(KIND=r8), INTENT(IN   ) :: fimax  ! global ! maximum fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: fi     (npoi)! global ! fractional snow cover
    REAL(KIND=r8), INTENT(INOUT) :: tsno   (npoi,nsnolay)! global ! temperature of snow layers (K)
    REAL(KIND=r8), INTENT(INOUT) :: hsno   (npoi,nsnolay)! global ! thickness of snow layers (m)

    !      INCLUDE 'comsoi.h'
    REAL(KIND=r8), INTENT(IN   ) :: sand    (npoi,nsoilay)! global ! percent sand of soil
    REAL(KIND=r8), INTENT(IN   ) :: clay    (npoi,nsoilay)! global ! percent clay of soil
    REAL(KIND=r8), INTENT(IN   ) :: poros   (npoi,nsoilay)! global ! porosity (mass of h2o per unit vol at sat / rhow)
    REAL(KIND=r8), INTENT(INOUT) :: wsoi    (npoi,nsoilay)! global ! fraction of soil pore space containing liquid water
    REAL(KIND=r8), INTENT(INOUT) :: wisoi   (npoi,nsoilay)! global ! fraction of soil pore space containing ice
    REAL(KIND=r8), INTENT(INOUT) :: consoi  (npoi,nsoilay)! local  ! thermal conductivity of each soil layer (W m-1 K-1)
    REAL(KIND=r8), INTENT(IN   ) :: zwpmax                   ! global ! assumed maximum fraction of soil surface 
    ! covered by puddles (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: wpud    (npoi)! global ! liquid content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: wipud   (npoi)! global ! ice content of puddles per soil area (kg m-2)
    REAL(KIND=r8), INTENT(IN   ) :: wpudmax                        ! normalization constant for puddles (kg m-2)
    REAL(KIND=r8), INTENT(INOUT) :: qglif   (npoi,4) ! local ! 1: fraction of soil evap (fvapg) from soil liquid
    ! 2: fraction of soil evap (fvapg) from soil ice
    ! 3: fraction of soil evap (fvapg) from puddle liquid
    ! 4: fraction of soil evap (fvapg) from puddle ice
    REAL(KIND=r8), INTENT(INOUT) :: tsoi    (npoi,nsoilay)! global        ! soil temperature for each layer (K)
    REAL(KIND=r8), INTENT(INOUT) :: hvasug  (npoi) ! local ! latent heat of vap/subl, for soil surface (J kg-1)
    REAL(KIND=r8), INTENT(INOUT) :: hvasui  (npoi)! local ! latent heat of vap/subl, for snow surface (J kg-1)
    REAL(KIND=r8), INTENT(IN   ) :: albsav  (npoi)! global ! saturated soil surface albedo (visible waveband)
    REAL(KIND=r8), INTENT(IN   ) :: albsan  (npoi)! global ! saturated soil surface albedo (near-ir waveband)
    REAL(KIND=r8), INTENT(INOUT) :: tg      (npoi)! global ! soil skin temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: ti      (npoi)! global ! snow skin temperature (K)
    REAL(KIND=r8), INTENT(IN   ) :: z0soi   (npoi)! global ! roughness length of soil surface (m)
    REAL(KIND=r8), INTENT(IN   ) :: swilt   (npoi,nsoilay)! global ! wilting soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(IN   ) :: sfield  (npoi,nsoilay)! global ! field capacity soil moisture value (fraction of pore space)
    REAL(KIND=r8), INTENT(INOUT) :: stressl (npoi,nsoilay)! local ! soil moisture stress factor for the lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stressu (npoi,nsoilay)! local ! soil moisture stress factor for the upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstl(npoi) ! local ! sum of stressl over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: stresstu(npoi)! local ! sum of stressu over all 6 soil layers (dimensionless)
    REAL(KIND=r8), INTENT(IN   ) :: csoi    (npoi,nsoilay)! global ! specific heat of soil, no pore spaces (J kg-1 deg-1)
    REAL(KIND=r8), INTENT(IN   ) :: rhosoi  (npoi,nsoilay)! global ! soil density (without pores, not bulk) (kg m-3)
    REAL(KIND=r8), INTENT(IN   ) :: hsoi    (nsoilay+1)   ! global ! soil layer thickness (m)
    REAL(KIND=r8), INTENT(IN   ) :: suction (npoi,nsoilay)! global ! saturated matric potential (m-h2o)
    REAL(KIND=r8), INTENT(IN   ) :: bex     (npoi,nsoilay)! global ! exponent "b" in soil water potential
    REAL(KIND=r8), INTENT(INOUT) :: upsoiu  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: upsoil  (npoi,nsoilay)! local  ! soil water uptake from transpiration (kg_h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: heatg   (npoi)        ! local        ! net heat flux into soil surface (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: heati   (npoi)        ! local        ! net heat flux into snow surface (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: hydraul (npoi,nsoilay)! global ! saturated hydraulic conductivity (m/s)
    REAL(KIND=r8), INTENT(INOUT) :: porosflo(npoi,nsoilay)! global ! porosity after reduction by ice content
    INTEGER, INTENT(IN   ) :: ibex    (npoi,nsoilay)! global ! nint(bex), used for cpu speed
    REAL(KIND=r8), INTENT(IN   ) :: bperm   ! global! lower b.c. for soil profile drainage 
    ! (0.0 = impermeable; 1.0 = fully permeable)
    REAL(KIND=r8), INTENT(INOUT) :: hflo    (npoi,nsoilay+1)  ! downward heat transport through soil layers (W m-2)


    !   INCLUDE 'comatm.h'
    REAL(KIND=r8), INTENT(IN   ) :: ta     (npoi)         ! global ! air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: asurd  (npoi,nband) ! local  ! direct albedo of surface system
    REAL(KIND=r8), INTENT(INOUT) :: asuri  (npoi,nband) ! local  ! diffuse albedo of surface system 
    REAL(KIND=r8), INTENT(IN   ) :: coszen (npoi)         ! global ! cosine of solar zenith angle
    REAL(KIND=r8), INTENT(IN   ) :: solad  (npoi,nband) ! global ! direct downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: solai  (npoi,nband) ! global ! diffuse downward solar flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: fira   (npoi)         ! global ! incoming ir flux (W m-2)
    REAL(KIND=r8), INTENT(IN   ) :: raina  (npoi)         ! global ! rainfall rate (mm/s or kg m-2 s-1)
    REAL(KIND=r8), INTENT(IN   ) :: qa     (npoi)         ! global ! specific humidity (kg_h2o/kg_air)
    REAL(KIND=r8), INTENT(IN   ) :: psurf  (npoi)         ! global ! surface pressure (Pa)
    REAL(KIND=r8), INTENT(IN   ) :: snowa  (npoi)         ! global ! snowfall rate (mm/s or kg m-2 s-1 of water)
    REAL(KIND=r8), INTENT(IN   ) :: ua     (npoi)         ! global ! wind speed (m s-1)
    REAL(KIND=r8), INTENT(IN   ) :: o2conc                 ! global ! o2 concentration (mol/mol)
    REAL(KIND=r8), INTENT(IN   ) :: co2conc                 ! global ! co2 concentration (mol/mol)

    !   INCLUDE 'com1d.h'
    REAL(KIND=r8), INTENT(INOUT) :: fwetu    (npoi)         ! local ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(INOUT) :: rliqu    (npoi)         ! local ! proportion of fwetu due to liquid
    REAL(KIND=r8), INTENT(INOUT) :: fwets    (npoi)         ! local ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(INOUT) :: rliqs    (npoi)         ! local ! proportion of fwets due to liquid
    REAL(KIND=r8), INTENT(INOUT) :: fwetl    (npoi)         ! local ! fraction of lower canopy stem & leaf area wetted by
    ! intercepted liquid and/or snow
    REAL(KIND=r8), INTENT(INOUT) :: rliql    (npoi)         ! local ! proportion of fwetl due to liquid
    REAL(KIND=r8), INTENT(INOUT) :: solu     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by upper 
    ! canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: sols     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by upper 
    ! canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: soll     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by lower 
    ! canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: solg     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by unit 
    ! snow-free soil (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: soli     (npoi)         ! local ! solar flux (direct + diffuse) absorbed by unit 
    ! snow surface (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: scalcoefl(npoi,4)   ! local ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(INOUT) :: scalcoefu(npoi,4)   ! local ! term needed in upper canopy scaling
    INTEGER, INTENT(INOUT) :: indsol   (npoi)         ! local ! index of current strip for points with positive coszen
    REAL(KIND=r8), INTENT(INOUT) :: albsod   (npoi)         ! local ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(INOUT) :: albsoi   (npoi)         ! local ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8), INTENT(INOUT) :: albsnd   (npoi)         ! local ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(INOUT) :: albsni   (npoi)         ! local ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8), INTENT(OUT  ) :: relod    (npoi)         ! local ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reloi    (npoi)         ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reupd    (npoi)         ! local ! upward direct radiation per unit incident direct 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: reupi    (npoi)         ! local ! upward diffuse radiation per unit incident diffuse 
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: ablod    (npoi)         ! local ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: abloi    (npoi)         ! local ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: flodd    (npoi)         ! local ! downward direct radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: dummy    (npoi)         ! local ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8), INTENT(INOUT) :: flodi    (npoi)         ! local ! downward diffuse radiation per unit incident direct
    ! radiation on lower canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: floii    (npoi)         ! local ! downward diffuse radiation per unit incident 
    ! diffuse radiation on lower canopy
    REAL(KIND=r8), INTENT(INOUT) :: terml    (npoi,7)         ! local ! term needed in lower canopy scaling
    REAL(KIND=r8), INTENT(INOUT) :: termu    (npoi,7)         ! local ! term needed in upper canopy scaling
    REAL(KIND=r8), INTENT(INOUT) :: abupd    (npoi)         ! local ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: abupi    (npoi)         ! local ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8), INTENT(INOUT) :: fupdd    (npoi)         ! local ! downward direct radiation per unit incident direct
    ! beam on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: fupdi    (npoi)         ! local ! downward diffuse radiation per unit icident direct
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: fupii    (npoi)         ! local ! downward diffuse radiation per unit incident diffuse
    ! radiation on upper canopy (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol2d    (npoi)         ! local ! direct downward radiation  out of upper canopy 
    ! per unit vegetated (upper) area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol2i    (npoi)         ! local ! diffuse downward radiation out of upper
    ! canopy per unit vegetated (upper) area(W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol3d    (npoi)         ! local ! direct downward radiation  out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: sol3i    (npoi)         ! local ! diffuse downward radiation out of upper
    ! canopy + gaps per unit grid cell area (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firb     (npoi)         ! local ! net upward ir radiation at reference
    ! atmospheric level za (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firs     (npoi)         ! local ! ir radiation absorbed by upper canopy stems (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firu     (npoi)         ! local ! ir raditaion absorbed by upper canopy leaves (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firl     (npoi)         ! local ! ir radiation absorbed by lower canopy leaves and stems (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firg     (npoi)         ! local ! ir radiation absorbed by soil/ice (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: firi     (npoi)         ! local ! ir radiation absorbed by snow (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: snowg    (npoi)         ! local ! snowfall rate at soil level (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: tsnowg   (npoi)         ! local ! snowfall temperature at soil level (K) 
    REAL(KIND=r8), INTENT(INOUT) :: tsnowl   (npoi)         ! local ! snowfall temperature below upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: pfluxl   (npoi)         ! local ! heat flux on lower canopy leaves & stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: raing    (npoi)         ! local ! rainfall rate at soil level (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: traing   (npoi)         ! local ! rainfall temperature at soil level (K)
    REAL(KIND=r8), INTENT(INOUT) :: trainl   (npoi)         ! local ! rainfall temperature below upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: snowl    (npoi)         ! local ! snowfall rate below upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: tsnowu   (npoi)         ! local ! snowfall temperature above upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: pfluxu   (npoi)         ! local ! heat flux on upper canopy leaves due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: rainu    (npoi)         ! local  ! rainfall rate above upper canopy (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: trainu   (npoi)         ! local  ! rainfall temperature above upper canopy (K)
    REAL(KIND=r8), INTENT(INOUT) :: snowu    (npoi)         ! local ! snowfall rate above upper canopy (kg h2o m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: pfluxs   (npoi)         ! local ! heat flux on upper canopy stems due to intercepted h2o (W m-2)
    REAL(KIND=r8), INTENT(INOUT) :: rainl    (npoi)         ! local ! rainfall rate below upper canopy (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: tfac                 ! local ! (ps/p) ** (rair/cair) for atmospheric level  (const)
    REAL(KIND=r8):: rhoa     (npoi)         ! local ! air density at za (allowing for h2o vapor) (kg m-3)
    REAL(KIND=r8), INTENT(INOUT) :: cp       (npoi)         ! local ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
    REAL(KIND=r8), INTENT(INOUT) :: za       (npoi)         ! local ! height above the surface of atmospheric forcing (m)
    REAL(KIND=r8), INTENT(INOUT) :: bdl      (npoi)         ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for
    ! laower canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: dil      (npoi)         ! local ! inverse of momentum diffusion coefficient within lower canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: z3       (npoi)         ! local ! effective top of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z4       (npoi)         ! local ! effective bottom of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z34      (npoi)         ! local ! effective middle of the lower canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: exphl    (npoi)         ! local ! exp(lamda/2*(z3-z4)) for lower canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(INOUT) :: expl     (npoi)         ! local ! exphl**2
    REAL(KIND=r8), INTENT(INOUT) :: displ    (npoi)         ! local ! zero-plane displacement height for lower canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: bdu      (npoi)         ! local ! aerodynamic coefficient ([(tau/rho)/u**2] for upper
    ! canopy (A31/A30 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT) :: diu      (npoi)         ! local ! inverse of momentum diffusion coefficient within upper canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: z1       (npoi)         ! local ! effective top of upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z2       (npoi)         ! local ! effective bottom of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: z12      (npoi)         ! local ! effective middle of the upper canopy (for momentum) (m)
    REAL(KIND=r8), INTENT(INOUT) :: exphu    (npoi)         ! local ! exp(lamda/2*(z3-z4)) for upper canopy (A30 Pollard & Thompson)
    REAL(KIND=r8), INTENT(INOUT) :: expu     (npoi)         ! local ! exphu**2
    REAL(KIND=r8), INTENT(INOUT) :: dispu    (npoi)         ! local ! zero-plane displacement height for upper canopy (m)
    REAL(KIND=r8), INTENT(INOUT) :: alogg    (npoi)         ! local ! log of soil roughness
    REAL(KIND=r8), INTENT(INOUT) :: alogi    (npoi)         ! local ! log of snow roughness
    REAL(KIND=r8), INTENT(INOUT) :: alogav   (npoi)         ! local ! average of alogi and alogg 
    REAL(KIND=r8), INTENT(INOUT) :: alog4    (npoi)         ! local ! log (max(z4, 1.1*z0sno, 1.1*z0soi)) 
    REAL(KIND=r8), INTENT(INOUT) :: alog3    (npoi)         ! local ! log (z3 - displ)
    REAL(KIND=r8), INTENT(INOUT) :: alog2    (npoi)         ! local ! log (z2 - displ)
    REAL(KIND=r8), INTENT(INOUT) :: alog1    (npoi)         ! local ! log (z1 - dispu) 
    REAL(KIND=r8), INTENT(INOUT) :: aloga    (npoi)         ! local ! log (za - dispu) 
    REAL(KIND=r8), INTENT(INOUT) :: u2       (npoi)         ! local ! wind speed at level z2 (m s-1)
    REAL(KIND=r8), INTENT(INOUT) :: alogu    (npoi)         ! local ! log (roughness length of upper canopy)
    REAL(KIND=r8), INTENT(INOUT) :: alogl    (npoi)         ! local ! log (roughness length of lower canopy)
    REAL(KIND=r8), INTENT(INOUT) :: richl    (npoi)         ! local ! richardson number for air above upper canopy (z3 to z2)
    REAL(KIND=r8), INTENT(INOUT) :: straml   (npoi)         ! local ! momentum correction factor for stratif between
    ! upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: strahl   (npoi)         ! local ! heat/vap correction factor for stratif between
    ! upper & lower canopy (z3 to z2) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: richu    (npoi)         ! local ! richardson number for air between upper & lower canopy (z1 to za)
    REAL(KIND=r8), INTENT(INOUT)  :: stramu   (npoi)         ! local ! momentum correction factor for stratif above
    ! upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: strahu   (npoi)         ! local ! heat/vap correction factor for stratif above
    ! upper canopy (z1 to za) (louis et al.)
    REAL(KIND=r8), INTENT(INOUT)  :: u1       (npoi)         ! local ! wind speed at level z1 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u12      (npoi)         ! local ! wind speed at level z12 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u3       (npoi)         ! local ! wind speed at level z3 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u34      (npoi)         ! local ! wind speed at level z34 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: u4       (npoi)         ! local ! wind speed at level z4 (m s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: cu       (npoi)         ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3) for
    ! upper air region (z12 --> za) (A35 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT)  :: cl       (npoi)         ! local ! air transfer coefficient (*rhoa) (m s-1 kg m-3)
    ! between the 2 canopies (z34 --> z12) (A36 Pollard & Thompson 1995)
    REAL(KIND=r8), INTENT(INOUT)  :: sg       (npoi)         ! local ! air-soil transfer coefficient
    REAL(KIND=r8), INTENT(INOUT)  :: si       (npoi)         ! local ! air-snow transfer coefficient
    REAL(KIND=r8), INTENT(INOUT)  :: fwetux   (npoi)         ! local ! fraction of upper canopy leaf area wetted if dew forms
    REAL(KIND=r8), INTENT(INOUT)  :: fwetsx   (npoi)         ! local ! fraction of upper canopy stem area wetted if dew forms
    REAL(KIND=r8), INTENT(INOUT)  :: fwetlx   (npoi)         ! local ! fraction of lower canopy leaf and stem area wetted if dew forms
    REAL(KIND=r8), INTENT(INOUT)  :: fsena    (npoi)         ! local ! downward sensible heat flux between za & z12 at za (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fseng    (npoi)         ! local ! upward sensible heat flux between soil surface & air at z34 (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fseni    (npoi)         ! local ! upward sensible heat flux between snow surface & air at z34 (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fsenu    (npoi)         ! local ! sensible heat flux from upper canopy leaves to air (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fsens    (npoi)         ! local ! sensible heat flux from upper canopy stems to air (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fsenl    (npoi)         ! local ! sensible heat flux from lower canopy to air (W m-2)
    REAL(KIND=r8), INTENT(INOUT)  :: fvapa    (npoi)         ! local ! downward h2o vapor flux between za & z12 at za (kg m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaput   (npoi)         ! local ! h2o vapor flux (transpiration from dry parts) 
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaps    (npoi)         ! local ! h2o vapor flux (evaporation from wet surface)
    ! between upper canopy stems and air at z12 (kg m-2 s-1 / SAI lower canopy / fu)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaplw   (npoi)         ! local ! h2o vapor flux (evaporation from wet surface) 
    ! between lower canopy leaves & stems and air at z34 (kg m-2 s-1/ LAI lower canopy/ fl)
    REAL(KIND=r8), INTENT(INOUT)  :: fvaplt   (npoi)         ! local ! h2o vapor flux (transpiration) 
    ! between lower canopy & air at z34 (kg m-2 s-1 / LAI lower canopy / fl)
    REAL(KIND=r8), INTENT(INOUT)  :: fvapg    (npoi)         ! local ! h2o vapor flux (evaporation) between soil & air 
    ! at z34 (kg m-2 s-1/bare ground fraction)
    REAL(KIND=r8), INTENT(INOUT)  :: fvapi    (npoi)         ! local ! h2o vapor flux (evaporation) between snow & air at z34 (kg m-2 s-1 / fi )
    REAL(KIND=r8), INTENT(INOUT)  :: fvapuw   (npoi)         ! local ! h2o vapor flux (evaporation from wet parts)
    ! between upper canopy leaves and air at z12 (kg m-2 s-1/ LAI upper canopy/ fu)
    REAL(KIND=r8), INTENT(INOUT)  :: td       (npoi)      ! global! daily average temperature (K)

    REAL(KIND=r8), INTENT(IN   )  :: vzero   (npoi)         ! global! a real array of zeros, of length npoi


    INTEGER, INTENT(IN   ) :: ndaypy              ! global! number of days per year
    REAL(KIND=r8), INTENT(OUT  ) :: nppdummy (npoi,npft)! local ! canopy NPP before accounting for stem and root respiration
    REAL(KIND=r8), INTENT(INOUT) :: tgpp     (npoi,npft)! local ! instantaneous GPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tgpptot  (npoi)         ! local ! instantaneous gpp (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tnpp     (npoi,npft)! local ! instantaneous NPP for each pft (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(IN   ) :: cbiow    (npoi,npft)! global! carbon in woody biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(IN   ) :: sapfrac  (npoi)         ! global! fraction of woody biomass that is in sapwood
    REAL(KIND=r8), INTENT(IN   ) :: cbior    (npoi,npft)! global! carbon in fine root biomass pool (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: tnpptot  (npoi)         ! local ! instantaneous npp (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: tco2root (npoi)         ! local ! instantaneous fine co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(OUT  ) :: tneetot  (npoi)         ! local ! instantaneous net ecosystem exchange of co2 per timestep (kg_C m-2/timestep)
    REAL(KIND=r8), INTENT(INOUT) :: tco2mic  (npoi)         ! local ! instantaneous microbial co2 flux from soil (mol-CO2 / m-2 / second)
    REAL(KIND=r8), INTENT(INOUT) :: a10td    (npoi)     ! global! 10-day average daily air temperature (K)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancub (npoi)     ! global! 10-day average canopy photosynthesis rate - broadleaf (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancuc (npoi)     ! global! 10-day average canopy photosynthesis rate - conifer (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancls (npoi)     ! global! 10-day average canopy photosynthesis rate - shrubs (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl3 (npoi)     ! global! 10-day average canopy photosynthesis rate - c3 grasses (mol_co2 m-2 s-1)
    REAL(KIND=r8), INTENT(INOUT) :: a10ancl4 (npoi)     ! global! 10-day average canopy photosynthesis rate - c4 grasses (mol_co2 m-2 s-1)

    INTEGER, INTENT(INOUT) :: ndtimes(npoi)             ! global! counter for daily average calculations
    REAL(KIND=r8), INTENT(INOUT) :: adrain    (npoi)! global! daily average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsnow    (npoi)! global! daily average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adaet     (npoi)! global! daily average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adtrunoff (npoi)! global! daily average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adsrunoff (npoi)! global! daily average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: addrainage(npoi)! global! daily average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: adrh      (npoi)! global! daily average rh (percent)
    REAL(KIND=r8), INTENT(INOUT) :: adsnod    (npoi)! global! daily average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: adsnof    (npoi)! global! daily average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoi    (npoi)! global! daily average soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoi    (npoi)! global! daily average soil temperature (c)
    REAL(KIND=r8), INTENT(INOUT) :: adwisoi   (npoi)! global! daily average soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtlaysoi (npoi)! global! daily average soil temperature (c) of top layer
    REAL(KIND=r8), INTENT(INOUT) :: adwlaysoi (npoi)! global! daily average soil moisture of top layer(fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adwsoic   (npoi)! global! daily average soil moisture using root profile weighting (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: adtsoic   (npoi)! global! daily average soil temperature (c) using profile weighting
    REAL(KIND=r8), INTENT(INOUT) :: adco2mic  (npoi)! global! daily accumulated co2 respiration from microbes (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2root (npoi)! global! daily accumulated co2 respiration from roots (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2soi  (npoi)! global! daily accumulated co2 respiration from soil(total) (kg_C m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: adco2ratio(npoi)! global! ratio of root to total co2 respiration
    REAL(KIND=r8), INTENT(INOUT) :: adnmintot (npoi)! global! daily accumulated net nitrogen mineralization (kg_N m-2 /day)
    REAL(KIND=r8), INTENT(INOUT) :: decompl   (npoi)! global! litter decomposition factor              (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: decomps   (npoi)! global! soil organic matter decomposition factor      (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: tnmin     (npoi)! global! instantaneous nitrogen mineralization (kg_N m-2/timestep)

    INTEGER, INTENT(IN   ) :: ndaypm    (12)  ! global! number of days per month


    INTEGER, INTENT(INOUT) :: nmtimes        (npoi)           ! global! counter for monthly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: amrain        (npoi)     ! global! monthly average rainfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsnow        (npoi)     ! global! monthly average snowfall rate (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amaet        (npoi)     ! global! monthly average aet (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtrunoff  (npoi)     ! global! monthly average total runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amsrunoff  (npoi)     ! global! monthly average surface runoff (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amdrainage (npoi)     ! global! monthly average drainage (mm/day)
    REAL(KIND=r8), INTENT(INOUT) :: amtemp        (npoi)     ! global! monthly average air temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amqa        (npoi)     ! global! monthly average specific humidity (kg-h2o/kg-air)
    REAL(KIND=r8), INTENT(INOUT) :: amsolar        (npoi)     ! global! monthly average incident solar radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirup        (npoi)     ! global! monthly average upward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amirdown   (npoi)     ! global! monthly average downward ir radiation (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amsens        (npoi)     ! global! monthly average sensible heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlatent   (npoi)     ! global! monthly average latent heat flux (W/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlaiu        (npoi)     ! global! monthly average lai for upper canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amlail        (npoi)     ! global! monthly average lai for lower canopy (m**2/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: amtsoi        (npoi)     ! global! monthly average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: amwsoi        (npoi)     ! global! monthly average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amwisoi        (npoi)     ! global! monthly average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amvwc        (npoi)     ! global! monthly average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amawc        (npoi)     ! global! monthly average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amsnod        (npoi)     ! global! monthly average snow depth (m)
    REAL(KIND=r8), INTENT(INOUT) :: amsnof        (npoi)     ! global! monthly average snow fraction (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: amnpp        (npoi,npft)! global! monthly total npp for each plant type (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amnpptot   (npoi)     ! local ! monthly total npp for ecosystem (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2mic   (npoi)     ! global! monthly total CO2 flux from microbial respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amco2root  (npoi)     ! global! monthly total CO2 flux from soil due to root
    ! respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2soi   (npoi)     ! local ! monthly total soil CO2 flux from microbial
    !          and root respiration (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(OUT  ) :: amco2ratio (npoi)     ! local ! monthly ratio of root to total co2 flux
    REAL(KIND=r8), INTENT(OUT  ) :: amneetot   (npoi)     ! local ! monthly total net ecosystem exchange of CO2 (kg-C/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amnmintot  (npoi)     ! global! monthly total N mineralization from microbes (kg-N/m**2/month)
    REAL(KIND=r8), INTENT(INOUT) :: amalbedo   (npoi)         
    REAL(KIND=r8), INTENT(INOUT) :: amtsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwsoil    (npoi, nsoilay) 
    REAL(KIND=r8), INTENT(INOUT) :: amwisoil   (npoi, nsoilay)
    REAL(KIND=r8), INTENT(INOUT) :: amts2      (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransu   (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: amtransl   (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: amsuvap    (npoi)     ! global
    REAL(KIND=r8), INTENT(INOUT) :: aminvap    (npoi)     ! global
    INTEGER, INTENT(INOUT) :: nytimes         (npoi)           ! global! counter for yearly average calculations
    REAL(KIND=r8), INTENT(INOUT) :: aysolar        (npoi)     ! global! annual average incident solar radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirup        (npoi)     ! global! annual average upward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayirdown   (npoi)     ! global! annual average downward ir radiation (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aysens        (npoi)     ! global! annual average sensible heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aylatent   (npoi)     ! global! annual average latent heat flux (w/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayprcp        (npoi)     ! global! annual average precipitation (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayaet        (npoi)     ! global! annual average aet (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrans        (npoi)     ! global! annual average transpiration (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aytrunoff  (npoi)     ! global! annual average total runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aysrunoff  (npoi)     ! global! annual average surface runoff (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydrainage (npoi)     ! global! annual average drainage (mm/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aydwtot        (npoi)     ! global! annual average soil+vegetation+snow water 
    ! recharge (mm/yr or kg_h2o/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aywsoi        (npoi)     ! global! annual average 1m soil moisture (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aywisoi        (npoi)     ! global! annual average 1m soil ice (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aytsoi        (npoi)     ! global! annual average 1m soil temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: ayvwc        (npoi)     ! global! annual average 1m volumetric water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: ayawc        (npoi)     ! global! annual average 1m plant-available water content (fraction)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstu (npoi)     ! global! annual average soil moisture stress 
    ! parameter for upper canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aystresstl(npoi)      ! global! annual average soil moisture stress 
    ! parameter for lower canopy (dimensionless)
    REAL(KIND=r8), INTENT(INOUT) :: aygpp     (npoi,npft) ! global! annual gross npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aygpptot  (npoi)      ! local ! annual total gpp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: aynpp     (npoi,npft) ! global! annual total npp for each plant type(kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: aynpptot  (npoi)      ! local ! annual total npp for ecosystem (kg-c/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2mic  (npoi)      ! global! annual total CO2 flux from microbial respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayco2root (npoi)      ! global! annual total CO2 flux from soil due to root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayco2soi  (npoi)      ! local ! annual total soil CO2 flux from microbial and 
    ! root respiration (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(OUT  ) :: ayneetot  (npoi)      ! local! annual total NEE for ecosystem (kg-C/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayrootbio (npoi)      ! global! annual average live root biomass (kg-C / m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynmintot (npoi)      ! global! annual total nitrogen mineralization (kg-N/m**2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ayalit    (npoi)      ! global! aboveground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayblit    (npoi)      ! global! belowground litter (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycsoi    (npoi)      ! global! total soil carbon (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aycmic    (npoi)      ! global! total soil carbon in microbial biomass (kg-c/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayanlit   (npoi)      ! global! aboveground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aybnlit   (npoi)      ! global! belowground litter nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: aynsoi    (npoi)      ! global! total soil nitrogen (kg-N/m**2)
    REAL(KIND=r8), INTENT(INOUT) :: ayalbedo  (npoi)  
    REAL(KIND=r8), INTENT(INOUT) :: totalit   (npoi)           ! global! total standing aboveground litter (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrlit   (npoi)           ! global! total root litter carbon belowground (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcsoi   (npoi)           ! global! total carbon in all soil pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totcmic   (npoi)           ! global! total carbon residing in microbial pools (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totanlit  (npoi)           ! global! total standing aboveground nitrogen in litter (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totrnlit  (npoi)      ! global! total root litter nitrogen belowground (kg_N m-2)
    REAL(KIND=r8), INTENT(INOUT) :: totnsoi   (npoi)      ! global! total nitrogen in soil (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totnmic   (npoi)      ! local! total nitrogen residing in microbial pool (kg_N m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totlit    (npoi)      ! local! total carbon in all litter pools (kg_C m-2)
    REAL(KIND=r8), INTENT(OUT  ) :: totfall   (npoi)      ! local! total litterfall and root turnover (kg_C m-2/year)
    REAL(KIND=r8), INTENT(OUT  ) :: totnlit   (npoi)      ! local! total nitrogen in all litter pools (kg_N m-2)

    REAL(KIND=r8), INTENT(INOUT) :: firefac   (npoi)     ! global! factor that respresents the annual average
    REAL(KIND=r8), INTENT(INOUT) :: wtot      (npoi)     ! global! total amount of water stored in snow, soil,
    ! puddels, and on vegetation (kg_h2o)
    ! fuel dryness of a grid cell, and hence characterizes the readiness to burn

    REAL(KIND=r8), INTENT(INOUT) :: storedn (npoi)        ! global ! total storage of N in soil profile (kg_N m-2) 
    REAL(KIND=r8), INTENT(INOUT) :: yrleach (npoi)        ! global ! annual total amount C leached from soil profile (kg_C m-2/yr)
    REAL(KIND=r8), INTENT(INOUT) :: ynleach (npoi)
    REAL(KIND=r8), INTENT(IN   ) :: falll   (npoi)     ! global ! annual leaf litter fall (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: fallr   (npoi)     ! global ! annual root litter input                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(IN   ) :: fallw   (npoi)     ! global! annual wood litter fall                    (kg_C m-2/year)
    REAL(KIND=r8), INTENT(INOUT) :: clitlm  (npoi)     ! global! carbon in leaf litter pool - metabolic       (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitls  (npoi)     ! global! carbon in leaf litter pool - structural      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrm  (npoi)     ! global! carbon in fine root litter pool - metabolic  (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrs  (npoi)     ! global! carbon in fine root litter pool - structural (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwm  (npoi)     ! global! carbon in woody litter pool - metabolic      (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitws  (npoi)     ! global! carbon in woody litter pool - structural     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislop(npoi)     ! global! carbon in soil - slow protected humus           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoislon(npoi)     ! global! carbon in soil - slow nonprotected humus     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: csoipas (npoi)     ! global! carbon in soil - passive humus                   (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitll  (npoi)     ! global! carbon in leaf litter pool - lignin           (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitrl  (npoi)     ! global! carbon in fine root litter pool - lignin     (kg_C m-2)
    REAL(KIND=r8), INTENT(INOUT) :: clitwl  (npoi)     ! global! carbon in woody litter pool - lignin           (kg_C m-2)


    REAL(KIND=r8), INTENT(IN   ) :: tc      (npoi)        ! global  ! coldest monthly temperature (C)
    REAL(KIND=r8), INTENT(INOUT) :: agddu   (npoi)        ! global  ! annual accumulated growing degree days for bud
    ! burst, upper canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: tempu   (npoi)        ! global  ! cold-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: agddl   (npoi)        ! global  ! annual accumulated growing degree days for bud burst,
    ! lower canopy (day-degrees)
    REAL(KIND=r8), INTENT(INOUT) :: templ   (npoi)        ! global  ! cold-phenology trigger for grasses/shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropu   (npoi)        ! global  ! drought-phenology trigger for trees (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropls  (npoi)        ! global  ! drought-phenology trigger for shrubs (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl4  (npoi)        ! global  ! drought-phenology trigger for c4 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(INOUT) :: dropl3  (npoi)        ! global  ! drought-phenology trigger for c3 grasses (non-dimensional)
    REAL(KIND=r8), INTENT(IN   ) :: plai    (npoi,npft)! global  ! total leaf area index of each plant functional type

    !
    ! Arguments (input)
    !
    INTEGER, INTENT(IN   ) :: iday         ! day number  (passed in)
    INTEGER, INTENT(IN   ) :: imonth         ! month number (passed in)
    INTEGER, INTENT(IN   ) :: iyear
    INTEGER, INTENT(IN   ) :: iyear0
    INTEGER, INTENT(IN   ) :: isimveg 
    INTEGER, INTENT(IN   ) :: spinmax 
    REAL(KIND=r8), INTENT(IN        ) :: ux  (npoi)
    REAL(KIND=r8), INTENT(IN        ) :: uy  (npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: taux(npoi)
    REAL(KIND=r8), INTENT(OUT  ) :: tauy(npoi)
    REAL(KIND=r8), INTENT(INOUT) :: ts2 (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: qs2 (npoi)
    REAL(KIND=r8), INTENT(INOUT) :: gdd0this(npoi)         
    REAL(KIND=r8), INTENT(INOUT) :: gdd5this(npoi) 
    REAL(KIND=r8), INTENT(OUT  ) :: bstar(npoi)
    !
    INTEGER :: ib  
    INTEGER :: nsol                  ! number of points in indsol
    INTEGER :: spin

    !
    ! Compute current time step zenith angle
    !
    !      call ibiszen (calday  ,loni    ,lati    ,coszen, kpti ,kptj)


    IF (mcsec == 0.0_r8) THEN
       ! 
       ! Calculates phenology once a day (beginning of the day)
       !

       CALL pheno(tc      , &! INTENT(IN        )
            agddu   , &! INTENT(INOUT) global
            tempu   , &! INTENT(INOUT) global
            agddl   , &! INTENT(INOUT) global
            templ   , &! INTENT(INOUT) global
            dropu   , &! INTENT(INOUT) global
            dropls  , &! INTENT(INOUT) global
            dropl4  , &! INTENT(INOUT) global
            dropl3  , &! INTENT(INOUT) global
            plai    , &! INTENT(IN        )
            frac    , &! INTENT(OUT  )
            lai     , &! INTENT(OUT  )
            fl      , &! INTENT(IN        )
            fu      , &! INTENT(IN        )
            zbot    , &! INTENT(OUT  )
            ztop    , &! INTENT(OUT  )
            a10td   , &! INTENT(IN        )
            a10ancub, &! INTENT(IN        )
            a10ancls, &! INTENT(IN        )
            a10ancl4, &! INTENT(IN        )
            a10ancl3, &! INTENT(IN        )
            td      , &! INTENT(IN        )
            npoi    , &! INTENT(IN        )
            npft    , &! INTENT(IN        )
            epsilon   )! INTENT(IN        )

       !
       IF (isimveg .EQ. 1) THEN
          !
          ! call soil biogeochemistry model
          !

          !
          ! Soil carbon acceleration model deleted in favor of spinmax
          ! specification at each restart (AAM - 3/14/02)
          !

          !          if (soicspin .eq. 1) then
          !
          !             if ((iyear - iyear0) .le.
          !    >          (spinfrac * (nspinsoil - iyear0 - eqyears))) then
          !                spinmax = int(spincons)
          !
          !             else if ((iyear - iyear0) .lt.
          !    >              (nspinsoil - iyear0 -  eqyears)) then
          !
          !                slope   = spincons / ((nspinsoil - iyear0 - eqyears) -
          !    >                    (spinfrac * (nspinsoil - iyear0 - eqyears)))
          !
          !                spinmax = int (spincons - (slope * ((iyear - iyear0) -
          !    >                (spinfrac * (nspinsoil - iyear0 - eqyears)))))
          !
          !                spinmax = max(spinmax,1)
          !
          !             else
          !
          !                spinmax = 1
          !
          !             endif            ! if (iyear - iyear0) ....
          !
          !          else 
          !
          !             spinmax = 1
          !
          !          endif               ! if (soicspin = 1)

          DO  spin = 1, spinmax
             CALL soilbgc (iyear    , &! INTENT(IN   )
                  iyear0   , &! INTENT(IN   )
                  imonth   , &! INTENT(IN   )
                  iday     , &! INTENT(IN   )
                  spin     , &! INTENT(IN   )
                  spinmax  , &! INTENT(IN   )
                  ayprcp   , &! INTENT(IN   )
                  falll    , &! INTENT(IN   )
                  fallr    , &! INTENT(IN   )
                  fallw    , &! INTENT(IN   )
                  clitlm   , &! INTENT(INOUT)
                  clitls   , &! INTENT(INOUT)
                  clitrm   , &! INTENT(INOUT)
                  clitrs   , &! INTENT(INOUT)
                  clitwm   , &! INTENT(INOUT)
                  clitws   , &! INTENT(INOUT)
                  csoislop , &! INTENT(INOUT)
                  csoislon , &! INTENT(INOUT)
                  csoipas  , &! INTENT(INOUT)
                  totcmic  , &! INTENT(INOUT)
                  clitll   , &! INTENT(INOUT)
                  clitrl   , &! INTENT(INOUT)
                  clitwl   , &! INTENT(INOUT)
                  decomps  , &! INTENT(IN   )
                  decompl  , &! INTENT(IN   )
                  tnmin    , &! INTENT(OUT  )
                  totnmic  , &! INTENT(OUT  )
                  totlit   , &! INTENT(OUT  )
                  totalit  , &! INTENT(OUT  )
                  totrlit  , &! INTENT(OUT  )
                  totcsoi  , &! INTENT(OUT  )
                  totfall  , &! INTENT(OUT  )
                  totnlit  , &! INTENT(OUT  )
                  totanlit , &! INTENT(OUT  )
                  totrnlit , &! INTENT(OUT  )
                  totnsoi  , &! INTENT(OUT  )
                  tco2mic  , &! INTENT(OUT  )
                  storedn  , &! INTENT(INOUT)
                  yrleach  , &! INTENT(INOUT)
                  ynleach  , &! INTENT(INOUT)
                  hsoi     , &! INTENT(IN   )
                  sand     , &! INTENT(IN   )
                  clay     , &! INTENT(IN   )
                  npoi     , &! INTENT(IN   )
                  nsoilay  , &! INTENT(IN   )
                  ndaypy     )! INTENT(IN   )

          END DO
          !
       END IF                  ! if (isimveg = 1)
       !
    END IF                    ! if (msec < dtime)

    !call the land surface model

    CALL lsxmain(ginvap       ,gsuvap       , & 
         gtrans       ,gtransu      , & 
         gtransl      ,grunof       , & 
         gdrain       ,gadjust      , & 
         a10scalparamu,a10daylightu , & 
         a10scalparaml,a10daylightl , & 
         vmax_pft     ,tau15        , & 
         kc15              ,ko15            , & 
         cimax        ,gammaub      , & 
         alpha3       ,theta3       , & 
         beta3        ,coefmub      , & 
         coefbub      ,gsubmin      , & 
         gammauc      ,coefmuc      , & 
         coefbuc      ,gsucmin      , & 
         gammals      ,coefmls      , & 
         coefbls      ,gslsmin      , & 
         gammal3      ,coefml3      , & 
         coefbl3      ,gsl3min      , & 
         gammal4      ,alpha4       , & 
         theta4       ,beta4        , & 
         coefml4      ,coefbl4      , & 
         gsl4min      ,wliqu        , & 
         wliqumax     ,wsnou        , & 
         wsnoumax     ,tu           , & 
         wliqs        ,wliqsmax     , & 
         wsnos        ,wsnosmax     , & 
         ts              ,wliql        , & 
         wliqlmax     ,wsnol        , & 
         wsnolmax     ,tl           , & 
         topparu      ,topparl      , & 
         fl              ,fu           , & 
         lai              ,sai          , & 
         rhoveg       ,tauveg       , & 
         orieh        ,oriev        , & 
         wliqmin      ,wsnomin      , & 
         t12              ,tdripu       , & 
         tblowu       ,tdrips       , & 
         tblows       ,t34            , & 
         tdripl       ,tblowl       , & 
         ztop              ,alaiml       , & 
         zbot              ,alaimu       , & 
         froot        ,q34            , & 
         q12              ,su            , & 
         cleaf        ,dleaf        , & 
         ss              ,cstem        , & 
         dstem        ,sl            , & 
         cgrass       ,ciub            , & 
         ciuc              ,exist        , & 
         csub              ,gsub            , & 
         csuc              ,gsuc            , & 
         agcub        ,agcuc        , & 
         ancub        ,ancuc        , & 
         totcondub    ,totconduc    , & 
         cils              ,cil3              , & 
         cil4              ,csls              , & 
         gsls              ,csl3              , & 
         gsl3              ,csl4              , & 
         gsl4              ,agcls        , & 
         agcl4        ,agcl3        , & 
         ancls        ,ancl4        , & 
         ancl3        ,totcondls    , & 
         totcondl3    ,totcondl4    , & 
         chu              ,chs              , & 
         chl              ,frac              , & 
         tlsub        , & 
         z0sno        , & 
         rhos              , & 
         consno       , & 
         hsnotop      , & 
         hsnomin      , & 
         fimin        , & 
         fimax        , & 
         fi              , & 
         tsno              , & 
         hsno              , & 
         sand              , & 
         clay              , & 
         poros        , & 
         wsoi              , & 
         wisoi        , & 
         consoi       , & 
         zwpmax       , & 
         wpud              , & 
         wipud        , & 
         wpudmax      , & 
         qglif        , & 
         tsoi              , & 
         hvasug       , & 
         hvasui       , & 
         albsav       , & 
         albsan       , & 
         tg              , & 
         ti              , & 
         z0soi        , & 
         swilt        , & 
         sfield       , & 
         stressl      , & 
         stressu      , & 
         stresstl     , & 
         stresstu     , & 
         csoi              , & 
         rhosoi       , & 
         hsoi              , & 
         suction      , & 
         bex              , & 
         upsoiu       , & 
         upsoil       , & 
         heatg        , & 
         heati        , & 
         hydraul      , & 
         porosflo     , & 
         ibex              , & 
         bperm        , & 
         hflo              , & 
         ta              , & 
         asurd        , & 
         asuri        , & 
         coszen       , & 
         solad        , & 
         solai        , & 
         fira              , & 
         raina        , & 
         qa              , & 
         psurf        , & 
         snowa        , & 
         ua              , & 
         o2conc       , & 
         co2conc      , & 
         fwetu        , & 
         rliqu        , & 
         fwets        , & 
         rliqs        , & 
         fwetl        , & 
         rliql        , & 
                                !nsol              , & 
         solu              , & 
         sols              , & 
         soll              , & 
         solg              , & 
         soli              , & 
         scalcoefl    , & 
         scalcoefu    , & 
         indsol       , & 
         albsod       , & 
         albsoi       , & 
         albsnd       , & 
         albsni       , & 
         relod        , & 
         reloi        , & 
         reupd        , & 
         reupi        , & 
         ablod        , & 
         abloi        , & 
         flodd        , & 
         dummy        , & 
         flodi        , & 
         floii        , & 
         terml        , & 
         termu        , & 
         abupd        , & 
         abupi        , & 
         fupdd        , & 
         fupdi        , & 
         fupii        , & 
         sol2d        , & 
         sol2i        , & 
         sol3d        , & 
         sol3i        , & 
         firb              , & 
         firs              , & 
         firu              , & 
         firl              , & 
         firg              , & 
         firi              , & 
         snowg        , & 
         tsnowg       , & 
         tsnowl       , & 
         pfluxl       , & 
         raing        , & 
         traing       , & 
         trainl       , & 
         snowl        , & 
         tsnowu       , & 
         pfluxu       , & 
         rainu        , & 
         trainu       , & 
         snowu        , & 
         pfluxs       , & 
         rainl        , & 
         tfac              , & 
         rhoa              , & 
         cp              , & 
         za              , & 
         bdl              , & 
         dil              , & 
         z3              , & 
         z4              , & 
         z34              , & 
         exphl        , & 
         expl              , & 
         displ        , & 
         bdu              , & 
         diu              , & 
         z1              , & 
         z2              , & 
         z12              , & 
         exphu        , & 
         expu              , & 
         dispu        , & 
         alogg        , & 
         alogi        , & 
         alogav       , & 
         alog4        , & 
         alog3        , & 
         alog2        , & 
         alog1        , & 
         aloga        , & 
         u2              , & 
         alogu        , & 
         alogl        , & 
         richl        , & 
         straml       , & 
         strahl       , & 
         richu        , & 
         stramu       , & 
         strahu       , & 
         u1              , & 
         u12              , & 
         u3              , & 
         u34              , & 
         u4              , & 
         cu              , & 
         cl              , & 
         sg              , & 
         si              , & 
         fwetux       , & 
         fwetsx       , & 
         fwetlx       , & 
         fsena        , & 
         fseng        , & 
         fseni        , & 
         fsenu        , & 
         fsens        , & 
         fsenl        , & 
         fvapa        , & 
         fvaput       , & 
         fvaps        , & 
         fvaplw       , & 
         fvaplt       , & 
         fvapg        , & 
         fvapi        , & 
         fvapuw       , & 
         npoi              , & 
         nband        , & 
         nsoilay      , & 
         nsnolay      , & 
         npft              , & 
         epsilon      , & 
         dtime        , & 
         stef              , & 
         vonk              , & 
         grav              , & 
         tmelt        , & 
         hfus              , & 
         hvap              , & 
         hsub              , & 
         ch2o              , & 
         cice              , & 
         cair              , & 
         cvap              , & 
         rair              , & 
         rvap              , & 
         cappa        , & 
         rhow              , & 
         vzero        , & 
         pi              , &
         ux              , &! INTENT(IN        ) !global
         uy              , &! INTENT(IN        ) !global
         taux              , &! INTENT(OUT        ) !local
         tauy              , &! INTENT(OUT        ) !local
         bstar        , &! INTENT(OUT        ) !local
         ts2              , &! INTENT(OUT        ) !local
         qs2          )  ! INTENT(OUT        ) !local

    !
    ! accumulate some variables every timestep
    !
    CALL  sumnow(a10td   , &! INTENT(INOUT) !global
         a10ancub, &! INTENT(INOUT) !global
         a10ancuc, &! INTENT(INOUT) !global
         a10ancls, &! INTENT(INOUT) !global
         a10ancl3, &! INTENT(INOUT) !global
         a10ancl4, &! INTENT(INOUT) !global
         nppdummy, &! INTENT(OUT  ) !local
         frac    , &! INTENT(IN   ) !global
         ancub   , &! INTENT(IN   ) !global
         lai           , &! INTENT(IN   ) !global
         fu           , &! INTENT(IN   ) !global
         ancuc   , &! INTENT(IN   ) !global
         ancls   , &! INTENT(IN   ) !global
         fl           , &! INTENT(IN   ) !global
         ancl4   , &! INTENT(IN   ) !global
         ancl3   , &! INTENT(IN   ) !global
         tgpp    , &! INTENT(OUT  ) !local
         agcub   , &! INTENT(IN   ) !global
         agcuc   , &! INTENT(IN   ) !global
         agcls   , &! INTENT(IN   ) !global
         agcl4   , &! INTENT(IN   ) !global
         agcl3   , &! INTENT(IN   ) !global
         tgpptot , &! INTENT(OUT  ) !local
         ts      , &! INTENT(IN   ) !global
         froot   , &! INTENT(IN   ) !global
         tnpp           , &! INTENT(OUT  ) !local
         cbiow   , &! INTENT(IN   ) !global
         sapfrac , &! INTENT(IN   ) !global
         cbior   , &! INTENT(IN   ) !global
         tnpptot , &! INTENT(OUT  ) !local
         tco2root, &! INTENT(OUT  ) !local
         tneetot , &! INTENT(OUT  ) !local
         tco2mic , &! INTENT(IN   ) !global
         tsoi           , &! INTENT(IN   ) !global
         fi           , &! INTENT(IN   ) !global
         td           , &! INTENT(IN   ) !global
         npoi    , &! INTENT(IN   ) !global
         nsoilay , &! INTENT(IN   ) !global
         npft           , &! INTENT(IN   ) !global
         ndaypy  , &! INTENT(IN   ) !global
         dtime     )! INTENT(IN   ) !global

    CALL sumday(raina     , &! INTENT(IN   )
         snowa     , &! INTENT(IN   )
         fvapa     , &! INTENT(IN   )
         grunof    , &! INTENT(IN   )
         gdrain    , &! INTENT(IN   )
         hsno      , &! INTENT(IN   )
         fi        , &! INTENT(IN   )
         hsoi      , &! INTENT(IN   )
         tsoi      , &! INTENT(IN   )
         wsoi      , &! INTENT(IN   )
         wisoi     , &! INTENT(IN   )
         ndtimes   , &! INTENT(INOUT) global
         adrain    , &! INTENT(INOUT) global
         adsnow    , &! INTENT(INOUT) global
         adaet     , &! INTENT(INOUT) global
         adtrunoff , &! INTENT(INOUT) global
         adsrunoff , &! INTENT(INOUT) global
         addrainage, &! INTENT(INOUT) global
         adrh      , &! INTENT(INOUT) global
         adsnod    , &! INTENT(INOUT) global
         adsnof    , &! INTENT(INOUT) global
         adwsoi    , &! INTENT(INOUT) global
         adtsoi    , &! INTENT(INOUT) global
         adwisoi   , &! INTENT(INOUT) global
         adtlaysoi , &! INTENT(INOUT) global
         adwlaysoi , &! INTENT(INOUT) global
         adwsoic   , &! INTENT(INOUT) global
         adtsoic   , &! INTENT(INOUT) global
         adco2mic  , &! INTENT(INOUT) global
         adco2root , &! INTENT(INOUT) global
         adco2soi  , &! INTENT(INOUT) global
         adco2ratio, &! INTENT(INOUT) global
         adnmintot , &! INTENT(INOUT) global
         froot     , &! INTENT(IN   )
         tco2mic   , &! INTENT(IN   )
         tco2root  , &! INTENT(IN   )
         decompl   , &! INTENT(INOUT) global
         decomps   , &! INTENT(INOUT) global
         tnmin     , &! INTENT(IN   )
         npoi      , &! INTENT(IN   )
         nsoilay   , &! INTENT(IN   )
         nsnolay   , &! INTENT(IN   )
         dtime     , &! INTENT(IN   )
         td            , &! INTENT(IN   )
         gdd0this  , &! INTENT(INOUT) global
         gdd5this  , &! INTENT(INOUT) global
         ts2       , &! INTENT(INOUT) global
         mcsec             ) ! INTENT(INOUT) global

    CALL summonth (dtime     , &! INTENT(IN   )!global
         mcsec     , &! INTENT(IN   )!global
         iday      , &! INTENT(IN   )!global
         imonth    , &! INTENT(IN   )!global
         nmtimes   , &! INTENT(INOUT)!global
         amrain    , &! INTENT(INOUT)!global
         amsnow    , &! INTENT(INOUT)!global
         amaet     , &! INTENT(INOUT)!global
         amtrunoff , &! INTENT(INOUT)!global
         amsrunoff , &! INTENT(INOUT)!global
         amdrainage, &! INTENT(INOUT)!global
         amtemp    , &! INTENT(INOUT)!global
         amqa            , &! INTENT(INOUT)!global
         amsolar   , &! INTENT(INOUT)!global
         amirup    , &! INTENT(INOUT)!global
         amirdown  , &! INTENT(INOUT)!global
         amsens    , &! INTENT(INOUT)!global
         amlatent  , &! INTENT(INOUT)!global
         amlaiu    , &! INTENT(INOUT)!global
         amlail    , &! INTENT(INOUT)!global
         amtsoi    , &! INTENT(INOUT)!global
         amwsoi    , &! INTENT(INOUT)!global
         amwisoi   , &! INTENT(INOUT)!global
         amvwc     , &! INTENT(INOUT)!global
         amawc     , &! INTENT(INOUT)!global
         amsnod    , &! INTENT(INOUT)!global
         amsnof    , &! INTENT(INOUT)!global
         amnpp            , &! INTENT(INOUT)!global
         amnpptot  , &! INTENT(OUT  )!local
         amco2mic  , &! INTENT(INOUT)!global
         amco2root , &! INTENT(INOUT)!global
         amco2soi  , &! INTENT(OUT  )!local
         amco2ratio, &! INTENT(OUT  )!local
         amneetot  , &! INTENT(OUT  )!local
         amnmintot , &! INTENT(INOUT)!global
         amts2     , &! INTENT(INOUT)!global
         amtransu  , &! INTENT(INOUT)!global
         amtransl  , &! INTENT(INOUT)!global
         amsuvap   , &! INTENT(INOUT)!global
         aminvap   , &! INTENT(INOUT)!global
         amalbedo  , &! INTENT(INOUT)!global
         amtsoil   , &! INTENT(INOUT)!global
         amwsoil   , &! INTENT(INOUT)!global
         amwisoil  , &! INTENT(INOUT)!global
         ts2       , &! INTENT(INOUT)!global
         fu        , &! INTENT(IN   )!global
         lai       , &! INTENT(IN   )!global
         fl            , &! INTENT(IN   )!global
         tnpp      , &! INTENT(IN   )!global
         tco2mic   , &! INTENT(IN   )!global
         tco2root  , &! INTENT(IN   )!global
         tnmin     , &! INTENT(IN   )!global
         hsoi      , &! INTENT(IN   )!global
         tsoi      , &! INTENT(IN   )!global
         wsoi      , &! INTENT(IN   )!global
         wisoi            , &! INTENT(IN   )!global
         poros     , &! INTENT(IN   )!global
         swilt            , &! INTENT(IN   )!global
         hsno      , &! INTENT(IN   )!global
         fi        , &! INTENT(IN   )!global
         grunof    , &! INTENT(IN   )!global
         gdrain    , &! INTENT(IN   )!global
         gtransu   , &! INTENT(IN   )!global
         gtransl   , &! INTENT(IN   )!global
         gsuvap    , &! INTENT(IN   )!global
         ginvap    , &! INTENT(IN   )!global
         asurd     , &! INTENT(IN   )!global
         asuri     , &! INTENT(IN   )!global
         fvapa     , &! INTENT(IN   )!global
         firb            , &! INTENT(IN   )!global
         fsena     , &! INTENT(IN   )!global
         raina     , &! INTENT(IN   )!global
         snowa     , &! INTENT(IN   )!global
         ta        , &! INTENT(IN   )!global
         qa        , &! INTENT(IN   )!global
         solad            , &! INTENT(IN   )!global
         solai     , &! INTENT(IN   )!global
         fira      , &! INTENT(IN   )!global
         npoi      , &! INTENT(IN   )!global
         nband            , &! INTENT(IN   )!global
         nsoilay   , &! INTENT(IN   )!global
         nsnolay   , &! INTENT(IN   )!global
         npft      , &! INTENT(IN   )!global
         ndaypm    , &! INTENT(IN   )!global
         hvap        )! INTENT(IN   )!global

    CALL sumyear(dtime     , &! INTENT(IN   )
         mcsec     , &! INTENT(IN   )
         iday      , &! INTENT(IN   )
         imonth    , &! INTENT(IN   )
         wliqu     , &! INTENT(IN   )
         wsnou     , &! INTENT(IN   )
         fu            , &! INTENT(IN   )
         lai            , &! INTENT(IN   )
         wliqs     , &! INTENT(IN   )
         wsnos     , &! INTENT(IN   )
         sai       , &! INTENT(IN   )
         wliql            , &! INTENT(IN   )
         wsnol     , &! INTENT(IN   )
         fl        , &! INTENT(IN   )
         tgpp      , &! INTENT(IN   )
         tnpp      , &! INTENT(IN   )
         firefac   , &! INTENT(INOUT) global
         tco2mic   , &! INTENT(IN   )
         tco2root  , &! INTENT(IN   )
         cbior     , &! INTENT(IN   )
         tnmin     , &! INTENT(IN   )
         totalit   , &! INTENT(IN   )
         totrlit   , &! INTENT(IN   )
         totcsoi   , &! INTENT(IN   )
         totcmic   , &! INTENT(IN   )
         totanlit  , &! INTENT(IN   )
         totrnlit  , &! INTENT(IN   )
         totnsoi   , &! INTENT(IN   )
         nytimes   , &! INTENT(INOUT) global
         aysolar   , &! INTENT(INOUT) global
         ayirup    , &! INTENT(INOUT) global
         ayirdown  , &! INTENT(INOUT) global
         aysens    , &! INTENT(INOUT) global
         aylatent  , &! INTENT(INOUT) global
         ayprcp    , &! INTENT(INOUT) global
         ayaet     , &! INTENT(INOUT) global
         aytrans   , &! INTENT(INOUT) global
         aytrunoff , &! INTENT(INOUT) global
         aysrunoff , &! INTENT(INOUT) global
         aydrainage, &! INTENT(INOUT) global
         aydwtot   , &! INTENT(INOUT) global
         aywsoi    , &! INTENT(INOUT) global
         aywisoi   , &! INTENT(INOUT) global
         aytsoi    , &! INTENT(INOUT) global
         ayvwc     , &! INTENT(INOUT) global
         ayawc     , &! INTENT(INOUT) global
         aystresstu, &! INTENT(INOUT) global
         aystresstl, &! INTENT(INOUT) global
         aygpp     , &! INTENT(INOUT) global
         aygpptot  , &! INTENT(OUT  ) local
         aynpp            , &! INTENT(INOUT) global
         aynpptot  , &! INTENT(OUT  ) local
         ayco2mic  , &! INTENT(INOUT) global
         ayco2root , &! INTENT(INOUT) global
         ayco2soi  , &! INTENT(OUT  ) global
         ayneetot  , &! INTENT(OUT  ) global
         ayrootbio , &! INTENT(INOUT) global
         aynmintot , &! INTENT(INOUT) global
         ayalit    , &! INTENT(INOUT) global
         ayblit    , &! INTENT(INOUT) global
         aycsoi    , &! INTENT(INOUT) global
         aycmic    , &! INTENT(INOUT) global
         ayanlit   , &! INTENT(INOUT) global
         aybnlit   , &! INTENT(INOUT) global
         aynsoi    , &! INTENT(INOUT) global
         ayalbedo  , &! INTENT(INOUT) global
         hsoi            , &! INTENT(IN   ) global
         wpud            , &! INTENT(IN   ) global
         wipud     , &! INTENT(IN   ) global
         poros     , &! INTENT(IN   ) global
         wsoi            , &! INTENT(IN   ) global
         wisoi     , &! INTENT(IN   ) global
         tsoi            , &! INTENT(IN   ) global
         swilt     , &! INTENT(IN   ) global
         stresstu  , &! INTENT(IN   ) global
         stresstl  , &! INTENT(IN   ) global
         fi            , &! INTENT(IN   ) global
         rhos      , &! INTENT(IN   ) global
         hsno            , &! INTENT(IN   ) global
         gtrans    , &! INTENT(IN   ) global
         grunof    , &! INTENT(IN   ) global
         gdrain    , &! INTENT(IN   ) global
         wtot            , &! INTENT(INOUT) global
         firb            , &! INTENT(IN   ) global
         fsena     , &! INTENT(IN   ) global
         fvapa     , &! INTENT(IN   ) global
         solad            , &! INTENT(IN   ) global
         solai     , &! INTENT(IN   ) global
         fira      , &! INTENT(IN   ) global
         raina     , &! INTENT(IN   ) global
         snowa            , &! INTENT(IN   ) global
         asurd     , &! INTENT(IN   ) global
         asuri     , &! INTENT(IN   ) global
         npoi            , &! INTENT(IN   ) global
         nband     , &! INTENT(IN   ) global
         nsoilay   , &! INTENT(IN   ) global
         nsnolay   , &! INTENT(IN   ) global
         npft      , &! INTENT(IN   ) global
         ndaypy    , &! INTENT(IN   ) global
         hvap      , &! INTENT(IN   ) global
         rhow        )! INTENT(IN   ) global
    IF (doalb) THEN
       !
       ! Compute next time step zenith angle
       !
       !         CALL ibiszen (calday1, loni, lati, coszen, kpti ,kptj)
       !
       ! Compute albedos (used in next step radiation computations)
       ! set up for solar calculations
       !
       !         CALL solset(loopi, kpti, kptj)
       !
       ! set up for solar calculations
       !
       CALL solset(npoi     , &! INTENT(IN   )
            nsol     , &! INTENT(OUT  )
            nband    , &! INTENT(IN   )
            solu     , &! INTENT(OUT  )
            sols     , &! INTENT(OUT  )
            soll     , &! INTENT(OUT  )
            solg     , &! INTENT(OUT  )
            soli     , &! INTENT(OUT  )
            scalcoefl, &! INTENT(OUT  )
            scalcoefu, &! INTENT(OUT  )
            indsol   , &! INTENT(OUT  )
            topparu  , &! INTENT(OUT  )
            topparl  , &! INTENT(OUT  )
            asurd    , &! INTENT(OUT  )
            asuri    , &! INTENT(OUT  )
            coszen     )! INTENT(IN   )  

       !
       ! solar calculations for each waveband
       !
       DO  ib = 1, nband
          !
          ! solsur sets surface albedos for soil and snow
          ! solalb performs the albedo calculations
          ! solarf uses the unit-incident-flux results from solalb
          ! to obtain absorbed fluxes sol[u,s,l,g,i] and 
          ! incident pars sunp[u,l]
          !
          CALL solsur (ib       , &! INTENT(IN   )
               tmelt    , &! INTENT(IN   )
               nsol     , &! INTENT(IN   )
               albsod   , &! INTENT(OUt  )
               albsoi   , &! INTENT(OUt  )
               albsnd   , &! INTENT(OUt  )
               albsni   , &! INTENT(OUt  )
               indsol   , &! INTENT(IN   )
               wsoi     , &! INTENT(IN   )
               wisoi    , &! INTENT(IN   )
               albsav   , &! INTENT(IN   )
               albsan   , &! INTENT(IN   )
               tsno     , &! INTENT(IN   )
               coszen   , &! INTENT(IN   )
               npoi     , &! INTENT(IN   )
               nsoilay  , &! INTENT(IN   )
               nsnolay    )! INTENT(IN   )

          CALL solalb (ib       , &! INTENT(IN   )
               relod    , &! INTENT(OUT  )
               reloi    , &! INTENT(OUT  )
               indsol   , &! INTENT(IN   )
               reupd    , &! INTENT(OUT  )
               reupi    , &! INTENT(OUT  )
               albsnd   , &! INTENT(IN   )
               albsni   , &! INTENT(IN   )
               albsod   , &! INTENT(IN   )
               albsoi   , &! INTENT(IN   )
               fl       , &! INTENT(IN   )
               fu       , &! INTENT(IN   )
               fi       , &! INTENT(IN   )
               asurd    , &! INTENT(INOUT)! local
               asuri    , &! INTENT(INOUT)! local
               npoi     , &! INTENT(IN   )
               nband    , &! INTENT(IN   )
               nsol     , &! INTENT(IN   )
               ablod    , &! INTENT(OUT  )
               abloi    , &! INTENT(OUT  )
               flodd    , &! INTENT(OUT  )
               dummy    , &! INTENT(OUT  )
               flodi    , &! INTENT(OUT  )
               floii    , &! INTENT(OUT  )
               coszen   , &! INTENT(IN   )
               terml    , &! INTENT(OUT  )
               termu    , &! INTENT(OUT  )
               lai      , &! INTENT(IN   )
               sai      , &! INTENT(IN   )
               abupd    , &! INTENT(OUT  )
               abupi    , &! INTENT(OUT  )
               fupdd    , &! INTENT(OUT  )
               fupdi    , &! INTENT(OUT  )
               fupii    , &! INTENT(OUT  )
               fwetl    , &! INTENT(IN   )
               rliql    , &! INTENT(IN   )
               rliqu    , &! INTENT(IN   )
               rliqs    , &! INTENT(IN   )
               fwetu    , &! INTENT(IN   )
               fwets    , &! INTENT(IN   )
               rhoveg   , &! INTENT(IN   )
               tauveg   , &! INTENT(IN   )
               orieh    , &! INTENT(IN   )
               oriev    , &! INTENT(IN   )
               tl       , &! INTENT(IN   )
               ts       , &! INTENT(IN   )
               tu       , &! INTENT(IN   )
               pi       , &! INTENT(IN   )
               tmelt    , &! INTENT(IN   )
               epsilon    )! INTENT(IN   )

          CALL solarf (ib       , & ! INTENT(IN        ) 
               nsol     , & ! INTENT(IN        ) 
               solu     , & ! INTENT(INOUT) !global
               indsol   , & ! INTENT(IN        ) 
               abupd    , & ! INTENT(IN        ) 
               abupi    , & ! INTENT(IN        ) 
               sols     , & ! INTENT(INOUT) !global
               sol2d    , & ! INTENT(OUT  ) 
               fupdd    , & ! INTENT(IN        ) 
               sol2i    , & ! INTENT(OUT  ) 
               fupii    , & ! INTENT(IN        ) 
               fupdi    , & ! INTENT(IN        ) 
               sol3d    , & ! INTENT(OUT  ) 
               sol3i    , & ! INTENT(OUT  ) 
               soll     , & ! INTENT(INOUT) !global
               ablod    , & ! INTENT(IN        ) 
               abloi    , & ! INTENT(IN        ) 
               flodd    , & ! INTENT(IN        ) 
               flodi    , & ! INTENT(IN        ) 
               floii    , & ! INTENT(IN        ) 
               solg     , & ! INTENT(INOUT) !global
               albsod   , & ! INTENT(IN        ) 
               albsoi   , & ! INTENT(IN        ) 
               soli     , & ! INTENT(INOUT) !global
               albsnd   , & ! INTENT(IN        ) 
               albsni   , & ! INTENT(IN        ) 
               scalcoefu, & ! INTENT(OUT  ) 
               termu    , & ! INTENT(IN        ) 
               scalcoefl, & ! INTENT(OUT  ) 
               terml    , & ! INTENT(IN        ) 
               lai      , & ! INTENT(IN        ) 
               sai      , & ! INTENT(IN        ) 
               fu       , & ! INTENT(IN        )
               fl       , & ! INTENT(IN        )
               topparu  , & ! INTENT(OUT  ) 
               topparl  , & ! INTENT(OUT  ) 
               solad    , & ! INTENT(IN        )
               solai    , & ! INTENT(IN        )
               npoi     , & ! INTENT(IN        ) 
               nband    , & ! INTENT(IN        ) 
               epsilon    ) ! INTENT(IN        ) 
          !
       END DO
       !
    END IF

  END SUBROUTINE Ibis







  SUBROUTINE seasfc_wgfs(tmtx  ,umtx  ,qmtx   , &
       slrad ,tsurf ,qsurf , &
       gu    ,gv    ,t_sib ,sh_sib,gps   ,tsea  ,dtc3x ,sinclt, &
       sigki ,delsig,sens  ,evap  ,umom  ,vmom  ,rmi   ,rhi   , &
       cond  ,stor  ,zorl  ,rnet  ,ncols ,kMax  ,Ustarm,z0    , &
       rho   ,qsfc  ,tsfc  ,mskant,bstar ,iMask ,HML   ,HUML  , &
       HVML  ,TSK   ,GSW   ,GLW)
    !
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle
    ! kpbl.......Number of layers pbl process is included( for u v,t )
    ! kqpbl......Number of layers pbl process is included( for q     )
    ! tmtx.......Temperature related matrix
    !            gmt(i,k,1)*d(gt(i,k-1))/dt+gmt(i,k,2)*d(gt(i,k))/dt=gmt(i,k,3)
    !            gmt(i,1,1)=0.
    !            gmt(*,*,1)...dimensionless
    !            gmt(*,*,2)...dimensionless
    !            gmt(*,*,3)...deg/sec
    ! umtx.......Wind related matrix
    !            gmu(i,k,1)*d(gu(i,k-1))/dt+gmu(i,k,2)*d(gu(i,k))/dt=gmu(i,k,3)
    !            gmu(i,k,1)*d(gv(i,k-1))/dt+gmu(i,k,2)*d(gv(i,k))/dt=gmu(i,k,4)
    !            gmu(i,1,1)=0.
    !            gmu(*,*,1)...dimensionless
    !            gmu(*,*,2)...dimensionless
    !            gmu(*,*,3)...m/sec**2
    !            gmu(*,*,4)...m/sec**2
    ! qmtx.......specific humidity related matrix
    !            gmq(i,k,1)*d(gq(i,k-1))/dt+gmq(i,k,2)*d(gq(i,k))/dt=gmq(i,k,3)
    !            gmq(i,1,1)=0.
    !            gmq(*,*,1)...dimensionless
    !            gmq(*,*,2)...dimensionless
    !            gmq(*,*,3)...kg/kg/sec
    ! slrad......radiation interpolation
    ! tsurff.....earth's surface temperature used for radiation
    !            for the first time step when ground temperature is not yet
    !            computed (this is done by subr.tsinit ),
    ! qsurf......qsurf(i)=0.622e0*EXP(21.65605e0 -5418.0e0 /tsurf(i))/gps(i)
    ! gu.........(zonal      velocity)*sin(colat)
    ! gv.........(meridional velocity)*sin(colat)
    ! gt.........Temperature
    ! gq.........Specific humidity
    ! gps........Surface pressure in mb
    ! tsea.......effective surface radiative temperature ( tgeff )
    ! dtc3x......time increment dt
    ! sinclt.....sinclt=SIN(colrad(latitu))
    ! sigki......sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !            sigma coordinate at middle of layer and akappa=gasr/cp
    ! delsig
    ! sens.......sensible heat flux
    ! evap.......latent heat flux  "evaporation"
    ! umom.......umom(i)=fmom*um(ncount),
    !            where .fmom  momentum flux      in n/m**2
    !            fmom= rhoair(ncount)*cu(ncount)*ustar(ncount)
    !            um  (ncount)=gu (i,1)/sinclt
    !            gu          = (zonal velocity)*sin(colat)
    ! vmom.......vmom(i)=rho(i)*gv(i)*rmi(i)
    !            rho  (i)=gps(i)/(gr100*gt(i))
    !            gr100 =gasr*0.01
    ! z0ice.......Roughness length of ice
    ! rmi.........rmi   (i)=cu(i)*ustar(i), where
    !             cu is friction  transfer coefficients
    !             ustar is surface friction velocity  (m/s)
    ! rhi.........rhi   (i)=ct(i)*ustar(i), where
    !             ct is heat transfer coefficients.
    !             ustar is surface friction velocity  (m/s)
    ! cond........cond(i)=gice*(tsurf(i)-tice) or
    !             cond(i)=(2.03/2.0)*(tsurf(i)-271.16)
    ! stor........stor(i)=hscap*c0(i)
    ! zorl........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !             zgrav =0.032 /grav
    ! rnet........rnet=-697.58*slrad(i)
    !             rnet(i)=rnet(i)-stefan*tsurf(i)**4
    ! cp..........specific heat of air           (j/kg/k)
    ! hl..........heat of evaporation of water     (j/kg)
    ! gasr........gas constant of dry air        (j/kg/k)
    ! grav........grav   gravity constant        (m/s**2)
    ! stefan......Stefan Boltzman constant
    !==========================================================================
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(INOUT) :: tmtx (ncols,3)
    REAL(KIND=r8),    INTENT(INOUT) :: umtx (ncols,4)
    REAL(KIND=r8),    INTENT(INOUT) :: qmtx (ncols,3)
    REAL(KIND=r8),    INTENT(IN   ) :: slrad(ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsurf(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: qsurf(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: gu   (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gv   (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: t_sib(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: sh_sib(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gps  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsea (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: dtc3x
    REAL(KIND=r8),    INTENT(IN   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: sigki
    REAL(KIND=r8),    INTENT(IN   ) :: delsig
    REAL(KIND=r8),    INTENT(INOUT) :: sens (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: evap (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: umom (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: vmom (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: rmi  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: rhi  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: cond (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: stor (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: zorl (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: rnet (ncols)
    REAL(KIND=r8),    INTENT(OUT  ) :: Ustarm  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: z0      (ncols)
    REAL(KIND=r8),    INTENT(OUT  ) :: rho   (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: qsfc (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsfc (ncols)
    INTEGER(KIND=i8), INTENT(IN   ) :: mskant(ncols)
    INTEGER(KIND=i8), INTENT(IN   ) :: iMask(ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: bstar(ncols)

    REAL(KIND=r8),    INTENT(INOUT) :: HML  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HUML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: HVML (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: TSK  (ncols)
    REAL(KIND=r8)    :: emisd(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: GSW(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: GLW(ncols)
    REAL(KIND=r8)    :: gt    (ncols) 
    REAL(KIND=r8)    :: gq    (ncols) 
    REAL(KIND=r8)    :: speedm  (ncols)
    REAL(KIND=r8)    :: ah    (ncols)
    REAL(KIND=r8)    :: al    (ncols)
    REAL(KIND=r8)    :: am    (ncols)
    REAL(KIND=r8)    :: cuni  (ncols)
    REAL(KIND=r8)    :: cui   (ncols)
    REAL(KIND=r8)    :: cu    (ncols)
    REAL(KIND=r8)    :: ctni  (ncols)
    REAL(KIND=r8)    :: cti   (ncols)
    REAL(KIND=r8)    :: ct    (ncols)
    REAL(KIND=r8)    :: um    (ncols)
    REAL(KIND=r8)    :: vm    (ncols)
    REAL(KIND=r8)    :: tha   (ncols)
    REAL(KIND=r8)    :: thm   (ncols)
    REAL(KIND=r8)    :: dzm   (ncols)
    REAL(KIND=r8)    :: thvgm (ncols)
    REAL(KIND=r8)    :: rib   (ncols)
    REAL(KIND=r8)    :: ustar (ncols)
    REAL(KIND=r8)    :: gtsav (ncols)
    REAL(KIND=r8)    :: gqsav (ncols)
    REAL(KIND=r8)    :: tmsav (ncols)
    REAL(KIND=r8)    :: qmsav (ncols)
    REAL(KIND=r8)    :: tssav (ncols)
    REAL(KIND=r8)    :: dqg0  (ncols)
    REAL(KIND=r8)    :: b00   (ncols)
    REAL(KIND=r8)    :: b03   (ncols)
    REAL(KIND=r8)    :: b04   (ncols)
    REAL(KIND=r8)    :: c0    (ncols)
    REAL(KIND=r8)    :: b30   (ncols)
    REAL(KIND=r8)    :: b33   (ncols)
    REAL(KIND=r8)    :: c3    (ncols)
    REAL(KIND=r8)    :: b40   (ncols)
    REAL(KIND=r8)    :: b44   (ncols)
    REAL(KIND=r8)    :: c4    (ncols)
    !LOGICAL :: jstneu
    !REAL(KIND=r8) :: u2 (ncols)

    INTEGER :: i
    INTEGER :: ncount
    REAL(KIND=r8)    :: gbyhl
    REAL(KIND=r8)    :: gbycp
    REAL(KIND=r8)    :: gr100
    REAL(KIND=r8)    :: gb100
    REAL(KIND=r8)    :: zgrav
    REAL(KIND=r8)    :: gice
    REAL(KIND=r8)    :: hscap
    REAL(KIND=r8)    :: sl1kap
    REAL(KIND=r8)    :: st4
    REAL(KIND=r8)    :: dti
    REAL(KIND=r8)    :: dtm
    REAL(KIND=r8)    :: dtmdt
    REAL(KIND=r8)    :: dqm
    REAL(KIND=r8)    :: dqmdt
    !*JPB REAL(KIND=r8), PARAMETER :: dd=0.05_r8
    REAL(KIND=r8), PARAMETER :: dd=3.0_r8 ! Total depth of the ice slab (m), Using ECMWF value
    REAL(KIND=r8), PARAMETER :: tice=271.16_r8
    REAL(KIND=r8), PARAMETER :: dice=2.0_r8
    REAL(KIND=r8), PARAMETER :: hice=2.03_r8
    REAL(KIND=r8), PARAMETER :: rhoice=920.0_r8 ! Mean ice density (kg/m3)
    REAL(KIND=r8), PARAMETER :: cice=2093.0_r8  ! Heat Capacity of Ice (J/Kg)
    REAL(KIND=r8) ::  U3D   (1:nCols)
    REAL(KIND=r8) ::  V3D   (1:nCols)
    REAL(KIND=r8) ::  T3D   (1:nCols)
    REAL(KIND=r8) ::  QV3D  (1:nCols)
    REAL(KIND=r8) ::  P3D   (1:nCols)
    REAL(KIND=r8) ::  PSFC  (1:nCols)
    REAL(KIND=r8) ::  CHS   (1:nCols)
    REAL(KIND=r8) ::  CHS2  (1:nCols)
    REAL(KIND=r8) ::  CQS2  (1:nCols)
    REAL(KIND=r8) ::  CPM   (1:nCols)
    !REAL(KIND=r8) ::  ZNT   (1:nCols)
    !REAL(KIND=r8) ::  UST   (1:nCols)
    REAL(KIND=r8) ::  PSIM  (1:nCols)
    REAL(KIND=r8) ::  PSIH  (1:nCols)
    REAL(KIND=r8) ::  XLAND (1:nCols)
    REAL(KIND=r8) ::  HFX   (1:nCols)
    REAL(KIND=r8) ::  QFX   (1:nCols)
    REAL(KIND=r8) ::  LH    (1:nCols)
    !REAL(KIND=r8) ::  TSK   (1:nCols)
    REAL(KIND=r8) ::  FLHC  (1:nCols)
    REAL(KIND=r8) ::  FLQC  (1:nCols)
    REAL(KIND=r8) ::  QGH   (1:nCols)
    REAL(KIND=r8) ::  QSFC_1  (1:nCols)
    REAL(KIND=r8) ::  U10   (1:nCols)
    REAL(KIND=r8) ::  V10   (1:nCols)
    REAL(KIND=r8) ::  GZ1OZ0(1:nCols)
    REAL(KIND=r8) ::  WSPD  (1:nCols)
    REAL(KIND=r8) ::  BR    (1:nCols)
    REAL(KIND=r8) :: CHS_SEA   (1:nCols)
    REAL(KIND=r8) :: CHS2_SEA  (1:nCols)
    REAL(KIND=r8) :: CPM_SEA   (1:nCols)
    REAL(KIND=r8) :: CQS2_SEA  (1:nCols)
    REAL(KIND=r8) :: FLHC_SEA  (1:nCols)
    REAL(KIND=r8) :: FLQC_SEA  (1:nCols)
    REAL(KIND=r8) :: HFX_SEA   (1:nCols)
    REAL(KIND=r8) :: LH_SEA    (1:nCols)
    REAL(KIND=r8) :: QFX_SEA   (1:nCols)
    REAL(KIND=r8) :: QGH_SEA   (1:nCols)
    REAL(KIND=r8) :: QSFC_SEA  (1:nCols)
    REAL(KIND=r8) :: UST_SEA   (1:nCols)
    REAL(KIND=r8) :: ZNT_SEA   (1:nCols)
    REAL(KIND=r8) ::  CM_SEA(1:nCols)
    REAL(KIND=r8) :: CH_SEA(1:nCols)
    REAL(KIND=r8) ::     WSPD_SEA(1:nCols)
    REAL(KIND=r8) :: SST       (1:nCols)
    REAL(KIND=r8) :: XICE      (1:nCols)
    REAL(KIND=r8) :: H0ML      (1:nCols)
    REAL(KIND=r8), PARAMETER ::  xice_threshold = 0.5_r8


    gr100 =rgas*0.01_r8
    gbycp =grav/(cp*delsig*100.0_r8 *sigki)
    gbyhl =grav/(hltm*delsig*100.0_r8 )
    gb100 =grav/(   delsig*100.0_r8 )
    zgrav =0.032_r8 /grav
    gice  =hice/dice ! 2.03_r8/2.0_r8
    hscap =rhoice*cice*dd/dtc3x
    sl1kap=sigki
    st4   =stefan*4.0_r8
    dti   =1.0_r8 /dtc3x
    gt=t_sib(:,1)
    gq=sh_sib(:,1)
   
    
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          rnet (i)=-697.58_r8*slrad(i)
          rho  (i)=gps(i)/(gr100*gt(i))
          ah   (i)=gbycp/gps(i)
          al   (i)=gbyhl/gps(i)
          dqg0 (i)=0.622_r8 *EXP(30.25353_r8 -5418.0_r8 /tsurf(i)) &
               /(tsurf(i)*tsurf(i)*gps(i))
          gtsav(i)=gt   (i)
          gqsav(i)=gq   (i)
          tssav(i)=tsurf(i)
          tmsav(i)=tmtx (i,3)
          qmsav(i)=qmtx (i,3)
       END IF
    END DO

    c0  =0.0_r8
    cond=0.0_r8
    stor=0.0_r8

    ncount=0
8000 CONTINUE
    ncount=ncount+1
    !
    !     the first call to vntlat just gets the neutral values of ustar
    !     and ventmf.
    !

    !    jstneu=.TRUE.
    !    CALL vntlt2 &
    !       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
    !       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
    !       thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,jstneu,u2  )
    !
    !   jstneu=.FALSE.
    !    CALL vntlt2 &
    !       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
    !       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
    !       thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,jstneu,u2  )
    !
    CALL vntlt1_wgfs ( &
         rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
         sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
         thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,kMax )

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt  (i)    =gtsav(i)
          gq  (i)    =gqsav(i)
          tsurf(i)   =tssav(i)
          tmtx(i,3)=tmsav(i)
          qmtx(i,3)=qmsav(i)
       END IF
    END DO
    
    DO i = 1, ncols   
       emisd(i) =0.98_r8   
       IF(mskant(i) == 1_i8)THEN
          U3D   (i) = gu (i,1)/sinclt(i)
          V3D   (i) = gv (i,1)/sinclt (i)        
          T3D   (i) = gt(i)
          QV3D  (i) = gq(i)
          P3D   (i) = gps(i)*100.0_r8 - gps(i)*delsig*100.0_r8
          PSFC  (i) = gps(i)*100.0_r8
          IF(iMask(i) >= 1_i8   )THEN ! land mask (1 for land, 2 for water , 15 ice)
            XLAND (i) =  2
          ELSE
            XLAND (i) =  2
         END IF
         IF(omlmodel)THEN
             H0ML  (i) = oml_hml0 - 13.5_r8*log(MAX(ABS(tsea(i))-tice+0.01_r8,1.0_r8)) 
             !H0ML  (i) = oml_hml0 + 2.0*sqrt(U3D(i)*U3D(i) + V3D(i)*U3D(i)) 
             IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
                TSK   (i) = ABS(tsurf(i))
             END IF
         ELSE 
            H0ML  (i) = 0.0_r8
            TSK   (i) = ABS(tsurf(i))
         END IF
         SST   (i) = ABS(tsea(i))!   REAL, INTENT(in ) ::, SST    (1:nCols)
         IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
            XICE  (i) = xice_threshold !   REAL, INTENT(in ) ::, XICE   (1:nCols)
         ELSE
            XICE  (i) = 0.0_r8
         END IF 
         z0    (i) = MAX(z0   (i),0.1_r8*z0ice)
      END IF
   END DO
    
    CALL sf_gfs_seaice_wrapper( &
       U3D       , &!   REAL, INTENT(IN   ) ::  U3D   (1:nCols) !-- U3D     3D u-velocity interpolated to theta points (m/s)
       V3D       , &!   REAL, INTENT(IN   ) ::  V3D   (1:nCols) !-- V3D     3D v-velocity interpolated to theta points (m/s)
       T3D       , &!   REAL, INTENT(IN   ) ::  T3D   (1:nCols) !-- T3D     temperature (K)
       QV3D      , &!   REAL, INTENT(IN   ) ::  QV3D  (1:nCols) !-- QV3D    3D water vapor mixing ratio (Kg/Kg)
       P3D       , &!   REAL, INTENT(IN   ) ::  P3D   (1:nCols) !-- P3D     3D pressure (Pa)
       PSFC      , &!   REAL, INTENT(IN   ) ::  PSFC  (1:nCols)                    surface pressure (Pa)
       CHS       , &!   REAL, INTENT(OUT  ) ::  CHS   (1:nCols)
       CHS2      , &!   REAL, INTENT(OUT  ) ::  CHS2  (1:nCols)
       CQS2      , &!   REAL, INTENT(OUT  ) ::  CQS2  (1:nCols)
       CPM       , &!   REAL, INTENT(OUT  ) ::  CPM   (1:nCols)
       z0        , &!   REAL, INTENT(INOUT) ::  ZNT   (1:nCols)
       ustar     , &!   REAL, INTENT(INOUT) ::  UST   (1:nCols)
       PSIM      , &!   REAL, INTENT(OUT  ) ::  PSIM  (1:nCols)
       PSIH      , &!   REAL, INTENT(OUT  ) ::  PSIH  (1:nCols)
       XLAND     , &!   REAL, INTENT(IN   ) ::  XLAND (1:nCols)
       HFX       , &!   REAL, INTENT(OUT  ) ::  HFX   (1:nCols)
       QFX       , &!   REAL, INTENT(OUT  ) ::  QFX   (1:nCols)
       LH        , &!   REAL, INTENT(OUT  ) ::  LH    (1:nCols)
       TSK       , &!   REAL, INTENT(IN   ) ::  TSK   (1:nCols)
       FLHC      , &!   REAL, INTENT(OUT  ) ::  FLHC  (1:nCols)
       FLQC      , &!   REAL, INTENT(OUT  ) ::  FLQC  (1:nCols)
       QGH       , &!   REAL, INTENT(OUT  ) ::  QGH    (1:nCols),
       QSFC_1    , &!   REAL, INTENT(OUT  ) ::  QSFC    (1:nCols),
       U10       , &!   REAL, INTENT(OUT  ) ::  U10     (1:nCols),
       V10       , &!   REAL, INTENT(OUT  ) ::  V10     (1:nCols),
       GZ1OZ0    , &!   REAL, INTENT(OUT  ) ::  GZ1OZ0  (1:nCols),
       WSPD      , &!   REAL, INTENT(OUT  ) ::  WSPD    (1:nCols),
       BR        , &!   REAL, INTENT(OUT  ) ::  BR      (1:nCols),
       CHS_SEA   , &!   REAL, INTENT(OUT) ::, CHS_SEA  (1:nCols)
       CHS2_SEA  , &!   REAL, INTENT(OUT) ::, CHS2_SEA  (1:nCols)
       CPM_SEA   , &!   REAL, INTENT(OUT) ::, CPM_SEA  (1:nCols)
       CQS2_SEA  , &!   REAL, INTENT(OUT) ::, CQS2_SEA  (1:nCols)
       FLHC_SEA  , &!   REAL, INTENT(OUT) ::, FLHC_SEA  (1:nCols)
       FLQC_SEA  , &!   REAL, INTENT(OUT) ::, FLQC_SEA  (1:nCols)
       HFX_SEA   , &!   REAL, INTENT(OUT) ::, HFX_SEA  (1:nCols)
       LH_SEA    , &!   REAL, INTENT(OUT) ::, LH_SEA  (1:nCols)
       QFX_SEA   , &!   REAL, INTENT(OUT) ::, QFX_SEA  (1:nCols)
       QGH_SEA   , &!   REAL, INTENT(OUT) ::, QGH_SEA  (1:nCols)
       QSFC_SEA  , &!   REAL, INTENT(OUT) ::, QSFC_SEA  (1:nCols)
       UST_SEA   , &!   REAL, INTENT(OUT) ::, UST_SEA  (1:nCols)
       ZNT_SEA   , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
       CM_SEA    , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
       CH_SEA    , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
       WSPD_SEA  , &!   REAL, INTENT(OUT) ::, ZNT_SEA  (1:nCols)
       SST       , &!   REAL, INTENT(in ) ::, SST    (1:nCols)
       XICE      , &!   REAL, INTENT(in ) ::, XICE   (1:nCols)
       mskant    , &
       delsig    , &
       nCols )
       IF(omlmodel)THEN
          CALL OCEANML( &
          ustar    , & !REAL,     INTENT(IN   ) :: UST  ( 1:nCols )
          U3D      , & !REAL,     INTENT(IN   ) :: U_PHY( 1:nCols )
          V3D      , & !REAL,     INTENT(IN   ) :: V_PHY( 1:nCols )
          mskant   , & !REAL,     INTENT(IN   ) :: XLAND( 1:nCols )
          HFX_SEA  , & !REAL,     INTENT(IN   ) :: HFX  ( 1:nCols )
          LH_SEA   , & !REAL,     INTENT(IN   ) :: LH   ( 1:nCols )
          ABS(tsea), & !REAL,     INTENT(IN   ) :: tsea ( 1:nCols )
          TSK      , & !REAL,     INTENT(INOUT) :: TSK  ( 1:nCols )
          HML      , & !REAL,     INTENT(INOUT) :: HML  ( 1:nCols )
          HUML     , & !REAL,     INTENT(INOUT) :: HUML ( 1:nCols )
          HVML     , & !REAL,     INTENT(INOUT) :: HVML ( 1:nCols )
          GSW      , & !REAL,     INTENT(IN   ) :: GSW  ( 1:nCols )
          GLW      , & !REAL,     INTENT(IN   ) :: GLW  ( 1:nCols )
          emisd    , & !REAL,     INTENT(IN   ) :: EMISS( 1:nCols )
          dtc3x    , & !REAL,     INTENT(IN   ) :: DT
          H0ML     , &
          nCols      ) !INTEGER,  INTENT(IN   ) :: nCols
       END IF

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             !
             ! Solution of sea ice
             !
             !
             !$$$$$$->PK
             !qsurf(i)=QSFC_SEA(i)
             !IF(omlmodel)THEN
             !   rhi(i)= CHS_SEA(i) 
             !   rmi  (i)=WSPD_SEA(i)*CM_SEA(i) 
             !END IF 
             !$$$$$$->PK
             !
             b00(i)=   hscap+cp*rho(i)*rhi(i) &
                  +hltm*rho(i)*rhi(i)*dqg0(i) &
                  +gice+st4*tsurf(i)**3
             b03(i)=        -cp*rho(i)*rhi(i)*sl1kap
             b04(i)=-hltm*rho(i)*rhi(i)
             !
             ! Right side of eq.41 section III.A 
             ! COLA Physics Description Manual
             !
             c0 (i)=rnet(i) -cp*rho(i)*rhi(i)*(tsurf(i)-sl1kap*gt(i)) &
                          -hltm*rho(i)*rhi(i)*(qsurf(i)        -       gq(i)) &
                  -gice*(tsurf(i)-tice)-stefan*tsurf(i)**4
             b30(i)=               -ah (i)*cp*rho(i)*rhi(i)
             b33(i)=tmtx(i,2)*dti-b30(i)*          sl1kap
             c3 (i)=tmtx(i,3)    -b30(i)*(tsurf(i)-sl1kap*gt(i))
             b40(i)=               -al(i)*hltm*rho(i)*rhi(i)* dqg0 (i)
             b44(i)=qmtx(i,2)*dti+al(i)*hltm*rho(i)*rhi(i)
             c4 (i)=qmtx(i,3)    + &
                  al(i)*hltm*rho(i)*rhi(i)*(qsurf(i)-gq(i))
             b00(i)=b00(i)-b30(i)*b03(i)/b33(i)-b40(i)*b04(i)/b44(i)
             c0 (i)=c0 (i)-c3 (i)*b03(i)/b33(i)-c4 (i)*b04(i)/b44(i)
             c0 (i)=c0 (i)/b00(i)
             tsurf(i)=tsurf(i)+c0(i)
             !IF(omlmodel)TSK  (i)=tsurf(i)+c0(i)
             tmtx(i,3)=(c3(i)-b30(i)*c0(i))/(b33(i)*dtc3x)
             qmtx(i,3)=(c4(i)-b40(i)*c0(i))/(b44(i)*dtc3x)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
             !
             ! Solution of sea water
             !
             !$$$$$$->PK
             rhi  (i) = CHS_SEA(i)
             rmi  (i) = WSPD_SEA(i)*(CM_SEA (i))
             !$$$$$$->PK
             zorl (i) = 100.0_r8    * zgrav*speedm(i)*rhi(i)
             !$$$$$$->PK
             sens (i) = HFX_SEA(i)
             evap (i) = LH_SEA(i)
             !pk  sens (i)= rho(i)*cp*(tsurf(i)-gt(i)*sigki(1))*rhi(i)
             !pk  evap (i)= rho(i)*hl*(qsurf(i)-gq(i))*rhi(i)
             !$$$$$$->PK
             tmtx(i,3)=(tmtx(i,3)+ah(i)*sens(i)) &
                  /(tmtx(i,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
             qmtx(i,3)=(qmtx(i,3)+al(i)*evap(i)) &
                  /(qmtx(i,2)+dtc3x*al(i)*rho(i)*hltm*rhi(i))
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt(i)=gt(i)+tmtx(i,3)*dtc3x
          gq(i)=gq(i)+qmtx(i,3)*dtc3x
       END IF
    END DO

    IF (ncount == 1) go to 8000

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN               
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             !
             ! Solution of sea ice
             !
             IF(omlmodel)THEN
                !rhi  (i)= CHS_SEA(i)
                !rmi  (i)= WSPD_SEA(i)*(CM_SEA(i))
                !sens (i) = HFX_SEA(i)
                !evap (i) = LH_SEA(i)
                sens(i)=rho(i)*cp*  (tsurf(i)   -  gt(i)*sigki   )*rhi(i)
                evap(i)=rho(i)*hltm*(qsurf(i)   -  gq(i)         )*rhi(i)
             ELSE
                sens(i)=rho(i)*cp*  (tsurf(i)   -  gt(i)*sigki   )*rhi(i)
                evap(i)=rho(i)*hltm*(qsurf(i)   -  gq(i)         )*rhi(i)
             END IF 
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
             !
             ! Solution of sea water
             !
             rhi  (i)= CHS_SEA(i)
             rmi  (i)= WSPD_SEA(i)*(CM_SEA(i))
             sens (i) = HFX_SEA(i)
             evap (i) = LH_SEA(i)
          END IF       
          dtmdt=(ah(i)*sens(i))/(tmtx(i,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i)+0.0001)
          dqmdt=(al(i)*evap(i))/(qmtx(i,2)+dtc3x*al(i)*rho(i)*hltm*rhi(i)+0.0001)
          dtm=dtmdt*dtc3x
          dqm=dqmdt*dtc3x
          tsfc   (i)=gt(i)+dtm
          qsfc   (i)=gq(i)+dqm

          gt  (i)=gtsav(i)
          gq  (i)=gqsav(i)
          rnet(i)=rnet(i)-stefan*tsurf(i)**4
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             cond(i)=gice*(tsurf(i)-tice)
             stor(i)=hscap*c0(i)
             rnet(i)=rnet(i)-stefan*tsurf(i)**3*4.0_r8 *c0(i)
             tsurf(i)=MIN(tsurf(i),tice)
             tsea (i)=-tsurf(i)
             !IF(omlmodel)TSK(i)  = tsurf(i)
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          !$$$$$$->PK
          !rmi  (i)= WSPD_SEA(i)*(CM_SEA(i) )
          !$$$$$$->PK
          umom(i)=rho(i)*gu(i,1)*rmi(i)
          vmom(i)=rho(i)*gv(i,1)*rmi(i)
          am  (i)=gb100/gps(i)
          umtx(i,3)=(umtx(i,3)-am(i)*umom(i)) &
               /(umtx(i,2)+dtc3x*am(i)*rho(i)*rmi(i))
          umtx(i,4)=(umtx(i,4)-am(i)*vmom(i)) &
               /(umtx(i,2)+dtc3x*am(i)*rho(i)*rmi(i))
          !
          !     set surface stress use of pseudo winds to true winds
          !     for output diagnostics
          !
          umom(i)=umom(i)/sinclt(i)
          vmom(i)=vmom(i)/sinclt(i)
          Ustarm(i) = SQRT(umom(i)**2 + vmom(i)**2)
          IF(Ustarm(i)==0.0_r8)Ustarm(i)=0.007_r8
          um  (i)=gu (i,1)/sinclt(i)
          vm  (i)=gv (i,1)/sinclt(i)
          speedm(i)=SQRT(um(i)**2 + vm(i)**2)
          speedm(i)=MAX(2.0_r8 , speedm(i))
          dzm   (i)=rbyg*gt(i)
          bstar(i)=cu(i)*grav*(ct(i)*(tsfc(i)-gt(i)-(grav/cp)*dzm(i))/gt(i))! + mapl_vireps*ct(i)*(qsfc(i)-gq(i)))
       END IF
    END DO
  END SUBROUTINE seasfc_wgfs




  SUBROUTINE seasfc_cola(tmtx  ,umtx  ,qmtx   , &
       slrad ,tsurf ,qsurf , &
       gu    ,gv    ,t_sib ,sh_sib ,gps  ,tsea  ,dtc3x ,sinclt, &
       sigki ,delsig,sens  ,evap  ,umom  ,vmom  ,rmi   ,rhi   , &
       cond  ,stor  ,zorl  ,rnet  ,ncols ,kMax  ,Ustarm,z0sea , &
       rho   ,qsfc  ,tsfc  ,mskant,bstar )
    !
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle
    ! kpbl.......Number of layers pbl process is included( for u v,t )
    ! kqpbl......Number of layers pbl process is included( for q     )
    ! tmtx.......Temperature related matrix
    !            gmt(i,k,1)*d(gt(i,k-1))/dt+gmt(i,k,2)*d(gt(i,k))/dt=gmt(i,k,3)
    !            gmt(i,1,1)=0.
    !            gmt(*,*,1)...dimensionless
    !            gmt(*,*,2)...dimensionless
    !            gmt(*,*,3)...deg/sec
    ! umtx.......Wind related matrix
    !            gmu(i,k,1)*d(gu(i,k-1))/dt+gmu(i,k,2)*d(gu(i,k))/dt=gmu(i,k,3)
    !            gmu(i,k,1)*d(gv(i,k-1))/dt+gmu(i,k,2)*d(gv(i,k))/dt=gmu(i,k,4)
    !            gmu(i,1,1)=0.
    !            gmu(*,*,1)...dimensionless
    !            gmu(*,*,2)...dimensionless
    !            gmu(*,*,3)...m/sec**2
    !            gmu(*,*,4)...m/sec**2
    ! qmtx.......specific humidity related matrix
    !            gmq(i,k,1)*d(gq(i,k-1))/dt+gmq(i,k,2)*d(gq(i,k))/dt=gmq(i,k,3)
    !            gmq(i,1,1)=0.
    !            gmq(*,*,1)...dimensionless
    !            gmq(*,*,2)...dimensionless
    !            gmq(*,*,3)...kg/kg/sec
    ! slrad......radiation interpolation
    ! tsurff.....earth's surface temperature used for radiation
    !            for the first time step when ground temperature is not yet
    !            computed (this is done by subr.tsinit ),
    ! qsurf......qsurf(i)=0.622e0*EXP(21.65605e0 -5418.0e0 /tsurf(i))/gps(i)
    ! gu.........(zonal      velocity)*sin(colat)
    ! gv.........(meridional velocity)*sin(colat)
    ! gt.........Temperature
    ! gq.........Specific humidity
    ! gps........Surface pressure in mb
    ! tsea.......effective surface radiative temperature ( tgeff )
    ! dtc3x......time increment dt
    ! sinclt.....sinclt=SIN(colrad(latitu))
    ! sigki......sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !            sigma coordinate at middle of layer and akappa=gasr/cp
    ! delsig
    ! sens.......sensible heat flux
    ! evap.......latent heat flux  "evaporation"
    ! umom.......umom(i)=fmom*um(ncount),
    !            where .fmom  momentum flux      in n/m**2
    !            fmom= rhoair(ncount)*cu(ncount)*ustar(ncount)
    !            um  (ncount)=gu (i,1)/sinclt
    !            gu          = (zonal velocity)*sin(colat)
    ! vmom.......vmom(i)=rho(i)*gv(i)*rmi(i)
    !            rho  (i)=gps(i)/(gr100*gt(i))
    !            gr100 =gasr*0.01
    ! z0ice.......Roughness length of ice
    ! rmi.........rmi   (i)=cu(i)*ustar(i), where
    !             cu is friction  transfer coefficients
    !             ustar is surface friction velocity  (m/s)
    ! rhi.........rhi   (i)=ct(i)*ustar(i), where
    !             ct is heat transfer coefficients.
    !             ustar is surface friction velocity  (m/s)
    ! cond........cond(i)=gice*(tsurf(i)-tice) or
    !             cond(i)=(2.03/2.0)*(tsurf(i)-271.16)
    ! stor........stor(i)=hscap*c0(i)
    ! zorl........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !             zgrav =0.032 /grav
    ! rnet........rnet=-697.58*slrad(i)
    !             rnet(i)=rnet(i)-stefan*tsurf(i)**4
    ! cp..........specific heat of air           (j/kg/k)
    ! hl..........heat of evaporation of water     (j/kg)
    ! gasr........gas constant of dry air        (j/kg/k)
    ! grav........grav   gravity constant        (m/s**2)
    ! stefan......Stefan Boltzman constant
    !==========================================================================
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(INOUT) :: tmtx (ncols,3)
    REAL(KIND=r8),    INTENT(INOUT) :: umtx (ncols,4)
    REAL(KIND=r8),    INTENT(INOUT) :: qmtx (ncols,3)
    REAL(KIND=r8),    INTENT(IN   ) :: slrad(ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsurf(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: qsurf(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: gu   (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gv   (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: t_sib(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: sh_sib(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gps  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsea (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: dtc3x
    REAL(KIND=r8),    INTENT(IN   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: sigki
    REAL(KIND=r8),    INTENT(IN   ) :: delsig
    REAL(KIND=r8),    INTENT(INOUT  ) :: sens (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: evap (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: umom (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: vmom (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: rmi  (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: rhi  (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: cond (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: stor (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: zorl (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: rnet (ncols)
    REAL(KIND=r8),    INTENT(OUT    ) :: Ustarm  (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: z0sea      (ncols)
    REAL(KIND=r8),    INTENT(OUT    ) :: rho   (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: qsfc (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: tsfc (ncols)
    INTEGER(KIND=i8), INTENT(IN     ) :: mskant(ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: bstar(ncols)


    REAL(KIND=r8)    :: z0    (ncols) 

    REAL(KIND=r8)    :: gt    (ncols) 
    REAL(KIND=r8)    :: gq    (ncols) 
    REAL(KIND=r8)    :: speedm  (ncols)
    REAL(KIND=r8)    :: ah    (ncols)
    REAL(KIND=r8)    :: al    (ncols)
    REAL(KIND=r8)    :: am    (ncols)
    REAL(KIND=r8)    :: cuni  (ncols)
    REAL(KIND=r8)    :: cui   (ncols)
    REAL(KIND=r8)    :: cu    (ncols)
    REAL(KIND=r8)    :: ctni  (ncols)
    REAL(KIND=r8)    :: cti   (ncols)
    REAL(KIND=r8)    :: ct    (ncols)
    REAL(KIND=r8)    :: um    (ncols)
    REAL(KIND=r8)    :: vm    (ncols)
    REAL(KIND=r8)    :: tha   (ncols)
    REAL(KIND=r8)    :: thm   (ncols)
    REAL(KIND=r8)    :: dzm   (ncols)
    REAL(KIND=r8)    :: thvgm (ncols)
    REAL(KIND=r8)    :: rib   (ncols)
    REAL(KIND=r8)    :: ustar (ncols)
    REAL(KIND=r8)    :: gtsav (ncols)
    REAL(KIND=r8)    :: gqsav (ncols)
    REAL(KIND=r8)    :: tmsav (ncols)
    REAL(KIND=r8)    :: qmsav (ncols)
    REAL(KIND=r8)    :: tssav (ncols)
    REAL(KIND=r8)    :: dqg0  (ncols)
    REAL(KIND=r8)    :: b00   (ncols)
    REAL(KIND=r8)    :: b03   (ncols)
    REAL(KIND=r8)    :: b04   (ncols)
    REAL(KIND=r8)    :: c0    (ncols)
    REAL(KIND=r8)    :: b30   (ncols)
    REAL(KIND=r8)    :: b33   (ncols)
    REAL(KIND=r8)    :: c3    (ncols)
    REAL(KIND=r8)    :: b40   (ncols)
    REAL(KIND=r8)    :: b44   (ncols)
    REAL(KIND=r8)    :: c4    (ncols)
    !LOGICAL :: jstneu
    !REAL(KIND=r8) :: u2 (ncols)

    INTEGER :: i
    INTEGER :: ncount
    REAL(KIND=r8)    :: gbyhl
    REAL(KIND=r8)    :: gbycp
    REAL(KIND=r8)    :: gr100
    REAL(KIND=r8)    :: gb100
    REAL(KIND=r8)    :: zgrav
    REAL(KIND=r8)    :: gice
    REAL(KIND=r8)    :: hscap
    REAL(KIND=r8)    :: sl1kap
    REAL(KIND=r8)    :: st4
    REAL(KIND=r8)    :: dti
    REAL(KIND=r8)    :: dtm
    REAL(KIND=r8)    :: dtmdt
    REAL(KIND=r8)    :: dqm
    REAL(KIND=r8)    :: dqmdt
    !*JPB REAL(KIND=r8), PARAMETER :: dd=0.05_r8
    REAL(KIND=r8), PARAMETER :: dd=3.0_r8 ! Total depth of the ice slab (m), Using ECMWF value
    REAL(KIND=r8), PARAMETER :: tice=271.16_r8
    REAL(KIND=r8), PARAMETER :: dice=2.0_r8
    REAL(KIND=r8), PARAMETER :: hice=2.03_r8
    REAL(KIND=r8), PARAMETER :: rhoice=920.0_r8 ! Mean ice density (kg/m3)
    REAL(KIND=r8), PARAMETER :: cice=2093.0_r8  ! Heat Capacity of Ice (J/Kg)


    sens =0.0_r8
    evap =0.0_r8
    umom =0.0_r8
    vmom =0.0_r8
    rmi  =0.0_r8
    rhi  =0.0_r8
    rnet =0.0_r8
    z0   =0.0_r8
    qsfc =0.0_r8
    tsfc =0.0_r8
    bstar=0.0_r8
    gt    =0.0_r8
    gq    =0.0_r8 
    speedm=0.0_r8
    ah    =0.0_r8
    al    =0.0_r8
    am    =0.0_r8
    cuni  =0.0_r8
    cui   =0.0_r8
    cu    =0.0_r8
    ctni  =0.0_r8
    cti   =0.0_r8
    ct    =0.0_r8
    um    =0.0_r8
    vm    =0.0_r8
    tha   =0.0_r8
    thm   =0.0_r8
    dzm   =0.0_r8
    thvgm =0.0_r8
    rib   =0.0_r8
    ustar =0.0_r8
    gtsav =0.0_r8
    gqsav =0.0_r8
    tmsav =0.0_r8
    qmsav =0.0_r8
    tssav =0.0_r8
    dqg0  =0.0_r8
    b00   =0.0_r8
    b03   =0.0_r8
    b04   =0.0_r8
    c0    =0.0_r8
    b30   =0.0_r8
    b33   =0.0_r8
    c3    =0.0_r8
    b40   =0.0_r8
    b44   =0.0_r8
    c4    =0.0_r8
    Ustarm=0.0_r8
    rho   =0.0_r8

    gbyhl =0.0_r8
    gbycp =0.0_r8
    gr100 =0.0_r8
    gb100 =0.0_r8
    zgrav =0.0_r8
    gice =0.0_r8
    hscap =0.0_r8
    sl1kap =0.0_r8
    st4 =0.0_r8
    dti =0.0_r8
    dtm =0.0_r8
    dtmdt =0.0_r8
    dqm =0.0_r8
    dqmdt =0.0_r8

    gr100 =rgas*0.01_r8
    gbycp =grav/(cp*delsig*100.0_r8 *sigki)
    gbyhl =grav/(hltm*delsig*100.0_r8 )
    gb100 =grav/(   delsig*100.0_r8 )
    zgrav =0.032_r8 /grav
    gice  =hice/dice ! 2.03_r8/2.0_r8
    hscap =rhoice*cice*dd/dtc3x
    sl1kap=sigki
    st4   =stefan*4.0_r8
    dti   =1.0_r8 /dtc3x
    gt=t_sib(:,1)
    gq=sh_sib(:,1)
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          rnet (i)=-697.58_r8*slrad(i)
          rho  (i)=gps(i)/(gr100*gt(i))
          ah   (i)=gbycp/gps(i)
          al   (i)=gbyhl/gps(i)
          dqg0 (i)=0.622_r8 *EXP(30.25353_r8 -5418.0_r8 /tsurf(i)) &
               /(tsurf(i)*tsurf(i)*gps(i))
          gtsav(i)=gt   (i)
          gqsav(i)=gq   (i)
          tssav(i)=tsurf(i)
          tmsav(i)=tmtx (i,3)
          qmsav(i)=qmtx (i,3)
       END IF
    END DO

    c0  =0.0_r8
    cond=0.0_r8
    stor=0.0_r8

    ncount=0
8000 CONTINUE
    ncount=ncount+1
    !
    !     the first call to vntlat just gets the neutral values of ustar
    !     and ventmf.
    !

    !    jstneu=.TRUE.
    !    CALL vntlt2 &
    !       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
    !       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
    !       thm   ,dzm   ,thvgm ,rib   ,z0sea ,zorl  ,ustar ,sinclt,mskant,jstneu,u2  )
    !
    !   jstneu=.FALSE.
    !    CALL vntlt2 &
    !       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
    !       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
    !       thm   ,dzm   ,thvgm ,rib   ,z0sea ,zorl  ,ustar ,sinclt,mskant,jstneu,u2  )
    !
    CALL vntlt1_cola ( &
         rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
         sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha    , &
         thm   ,dzm   ,thvgm ,rib   ,z0sea ,zorl  ,ustar ,sinclt,mskant ,kMax)

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt  (i)    =gtsav(i)
          gq  (i)    =gqsav(i)
          tsurf(i)   =tssav(i)
          tmtx(i,3)=tmsav(i)
          qmtx(i,3)=qmsav(i)
       END IF
    END DO
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             !
             ! Solution of sea ice
             !
             !
             b00(i)=   hscap+cp*rho(i)*rhi(i) &
                  +hltm*rho(i)*rhi(i)*dqg0(i) &
                  +gice+st4*tsurf(i)**3
             b03(i)=        -cp*rho(i)*rhi(i)*sl1kap
             b04(i)=-hltm*rho(i)*rhi(i)
             !
             ! Right side of eq.41 section III.A 
             ! COLA Physics Description Manual
             !
             c0 (i)=rnet(i) -cp*rho(i)*rhi(i)*(tsurf(i)-sl1kap*gt(i)) &
                  -hltm*rho(i)*rhi(i)*(qsurf(i)-       gq(i)) &
                  -gice*(tsurf(i)-tice)-stefan*tsurf(i)**4
             b30(i)=               -ah (i)*cp*rho(i)*rhi(i)
             b33(i)=tmtx(i,2)*dti-b30(i)*          sl1kap
             c3 (i)=tmtx(i,3)    -b30(i)*(tsurf(i)-sl1kap*gt(i))
             b40(i)=               -al(i)*hltm*rho(i)*rhi(i)* dqg0 (i)
             b44(i)=qmtx(i,2)*dti+al(i)*hltm*rho(i)*rhi(i)
             c4 (i)=qmtx(i,3)    + &
                  al(i)*hltm*rho(i)*rhi(i)*(qsurf(i)-gq(i))
             b00(i)=b00(i)-b30(i)*b03(i)/b33(i)-b40(i)*b04(i)/b44(i)
             c0 (i)=c0 (i)-c3 (i)*b03(i)/b33(i)-c4 (i)*b04(i)/b44(i)
             c0 (i)=c0 (i)/b00(i)
             tsurf(i)=tsurf(i)+c0(i)
             tmtx(i,3)=(c3(i)-b30(i)*c0(i))/(b33(i)*dtc3x)
             qmtx(i,3)=(c4(i)-b40(i)*c0(i))/(b44(i)*dtc3x)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
             !
             ! Solution of sea water
             !
             zorl (i)= 100.0_r8 *zgrav*speedm(i)*rhi(i)
             sens (i)= rho(i)*cp*(tsurf(i)-gt(i)*sigki)*rhi(i)
             evap (i)= rho(i)*hltm*(qsurf(i)-gq(i))*rhi(i)
             tmtx(i,3)=(tmtx(i,3)+ah(i)*sens(i)) &
                  /(tmtx(i,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
             qmtx(i,3)=(qmtx(i,3)+al(i)*evap(i)) &
                  /(qmtx(i,2)+dtc3x*al(i)*rho(i)*hltm*rhi(i))
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt(i)=gt(i)+tmtx(i,3)*dtc3x
          gq(i)=gq(i)+qmtx(i,3)*dtc3x
       END IF
    END DO

    IF (ncount == 1) go to 8000

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN          
          sens(i)=rho(i)*cp*(tsurf(i)-gt(i)*sigki)*rhi(i)
          evap(i)=rho(i)*hltm*(qsurf(i)-gq(i)         )*rhi(i)
          dtmdt=(ah(i)*sens(i))/(tmtx(i,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i)+0.0001)
          dqmdt=(al(i)*evap(i))/(qmtx(i,2)+dtc3x*al(i)*rho(i)*hltm*rhi(i)+0.0001)
          dtm=dtmdt*dtc3x
          dqm=dqmdt*dtc3x
          tsfc   (i)=gt(i)+dtm
          qsfc   (i)=gq(i)+dqm

          gt  (i)=gtsav(i)
          gq  (i)=gqsav(i)
          rnet(i)=rnet(i)-stefan*tsurf(i)**4
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             cond(i)=gice*(tsurf(i)-tice)
             stor(i)=hscap*c0(i)
             rnet(i)=rnet(i)-stefan*tsurf(i)**3*4.0_r8 *c0(i)
             tsurf(i)=MIN(tsurf(i),tice)
             tsea (i)=-   tsurf(i)
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          umom(i)=rho(i)*gu(i,1)*rmi(i)
          vmom(i)=rho(i)*gv(i,1)*rmi(i)
          am  (i)=gb100/gps(i)
          umtx(i,3)=(umtx(i,3)-am(i)*umom(i)) &
               /(umtx(i,2)+dtc3x*am(i)*rho(i)*rmi(i))
          umtx(i,4)=(umtx(i,4)-am(i)*vmom(i)) &
               /(umtx(i,2)+dtc3x*am(i)*rho(i)*rmi(i))
          !
          !     set surface stress use of pseudo winds to true winds
          !     for output diagnostics
          !
          umom(i)=umom(i)/sinclt(i)
          vmom(i)=vmom(i)/sinclt(i)
          Ustarm(i) = SQRT(umom(i)**2 + vmom(i)**2)
          IF(Ustarm(i)==0.0_r8)Ustarm(i)=0.007_r8
          um  (i)=gu (i,1)/sinclt(i)
          vm  (i)=gv (i,1)/sinclt(i)
          speedm(i)=SQRT(um(i)**2 + vm(i)**2)
          speedm(i)=MAX(2.0_r8 , speedm(i))
          dzm   (i)=rbyg*gt(i)
          bstar(i)=cu(i)*grav*(ct(i)*(tsfc(i)-gt(i)-(grav/cp)*dzm(i))/gt(i))! + mapl_vireps*ct(i)*(qsfc(i)-gq(i)))


       END IF
    END DO
  END SUBROUTINE seasfc_cola

  SUBROUTINE seasfc_ukme(tmtx  ,umtx  ,qmtx   , &
       slrad ,tsurf ,qsurf , &
       gu    ,gv    ,t_sib ,sh_sib,gps   ,tsea  ,dtc3x ,sinclt   , &
       sigki ,delsig,sens  ,evap  ,umom  ,vmom  ,rmi   ,rhi      , &
       cond  ,stor  ,zorl  ,rnet  ,ncols ,kMax  ,Ustarm,z0sea    , &
       rho   ,qsfc  ,tsfc  ,mskant,bstar ,iMask ,cldtot,ySwSfcNet, &
       LwSfcNet,pblh,QCF   ,QCL)
    !
    !==========================================================================
    ! ncols......Number of grid points on a gaussian latitude circle
    ! kpbl.......Number of layers pbl process is included( for u v,t )
    ! kqpbl......Number of layers pbl process is included( for q     )
    ! tmtx.......Temperature related matrix
    !            gmt(i,k,1)*d(gt(i,k-1))/dt+gmt(i,k,2)*d(gt(i,k))/dt=gmt(i,k,3)
    !            gmt(i,1,1)=0.
    !            gmt(*,*,1)...dimensionless
    !            gmt(*,*,2)...dimensionless
    !            gmt(*,*,3)...deg/sec
    ! umtx.......Wind related matrix
    !            gmu(i,k,1)*d(gu(i,k-1))/dt+gmu(i,k,2)*d(gu(i,k))/dt=gmu(i,k,3)
    !            gmu(i,k,1)*d(gv(i,k-1))/dt+gmu(i,k,2)*d(gv(i,k))/dt=gmu(i,k,4)
    !            gmu(i,1,1)=0.
    !            gmu(*,*,1)...dimensionless
    !            gmu(*,*,2)...dimensionless
    !            gmu(*,*,3)...m/sec**2
    !            gmu(*,*,4)...m/sec**2
    ! qmtx.......specific humidity related matrix
    !            gmq(i,k,1)*d(gq(i,k-1))/dt+gmq(i,k,2)*d(gq(i,k))/dt=gmq(i,k,3)
    !            gmq(i,1,1)=0.
    !            gmq(*,*,1)...dimensionless
    !            gmq(*,*,2)...dimensionless
    !            gmq(*,*,3)...kg/kg/sec
    ! slrad......radiation interpolation
    ! tsurff.....earth's surface temperature used for radiation
    !            for the first time step when ground temperature is not yet
    !            computed (this is done by subr.tsinit ),
    ! qsurf......qsurf(i)=0.622e0*EXP(21.65605e0 -5418.0e0 /tsurf(i))/gps(i)
    ! gu.........(zonal      velocity)*sin(colat)
    ! gv.........(meridional velocity)*sin(colat)
    ! gt.........Temperature
    ! gq.........Specific humidity
    ! gps........Surface pressure in mb
    ! tsea.......effective surface radiative temperature ( tgeff )
    ! dtc3x......time increment dt
    ! sinclt.....sinclt=SIN(colrad(latitu))
    ! sigki......sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !            sigma coordinate at middle of layer and akappa=gasr/cp
    ! delsig
    ! sens.......sensible heat flux
    ! evap.......latent heat flux  "evaporation"
    ! umom.......umom(i)=fmom*um(ncount),
    !            where .fmom  momentum flux      in n/m**2
    !            fmom= rhoair(ncount)*cu(ncount)*ustar(ncount)
    !            um  (ncount)=gu (i,1)/sinclt
    !            gu          = (zonal velocity)*sin(colat)
    ! vmom.......vmom(i)=rho(i)*gv(i)*rmi(i)
    !            rho  (i)=gps(i)/(gr100*gt(i))
    !            gr100 =gasr*0.01
    ! z0ice.......Roughness length of ice
    ! rmi.........rmi   (i)=cu(i)*ustar(i), where
    !             cu is friction  transfer coefficients
    !             ustar is surface friction velocity  (m/s)
    ! rhi.........rhi   (i)=ct(i)*ustar(i), where
    !             ct is heat transfer coefficients.
    !             ustar is surface friction velocity  (m/s)
    ! cond........cond(i)=gice*(tsurf(i)-tice) or
    !             cond(i)=(2.03/2.0)*(tsurf(i)-271.16)
    ! stor........stor(i)=hscap*c0(i)
    ! zorl........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !             zgrav =0.032 /grav
    ! rnet........rnet=-697.58*slrad(i)
    !             rnet(i)=rnet(i)-stefan*tsurf(i)**4
    ! cp..........specific heat of air           (j/kg/k)
    ! hl..........heat of evaporation of water     (j/kg)
    ! gasr........gas constant of dry air        (j/kg/k)
    ! grav........grav   gravity constant        (m/s**2)
    ! stefan......Stefan Boltzman constant
    !==========================================================================
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(INOUT) :: tmtx (ncols,3)
    REAL(KIND=r8),    INTENT(INOUT) :: umtx (ncols,4)
    REAL(KIND=r8),    INTENT(INOUT) :: qmtx (ncols,3)
    REAL(KIND=r8),    INTENT(IN   ) :: slrad(ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsurf(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: qsurf(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: gu   (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gv   (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: t_sib(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: sh_sib(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: gps  (ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: tsea (ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: dtc3x
    REAL(KIND=r8),    INTENT(IN   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: sigki
    REAL(KIND=r8),    INTENT(IN   ) :: delsig
    REAL(KIND=r8),    INTENT(INOUT  ) :: sens (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: evap (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: umom (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: vmom (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: rmi  (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: rhi  (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: cond (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: stor (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: zorl (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: rnet (ncols)
    REAL(KIND=r8),    INTENT(OUT    ) :: Ustarm  (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: z0sea      (ncols)
    REAL(KIND=r8),    INTENT(OUT    ) :: rho   (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: qsfc (ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: tsfc (ncols)
    INTEGER(KIND=i8), INTENT(IN     ) :: mskant(ncols)
    INTEGER(KIND=i8), INTENT(IN     ) :: iMask(ncols)
    REAL(KIND=r8),    INTENT(INOUT  ) :: bstar(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: cldtot (ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: ySwSfcNet(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: LwSfcNet(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: pblh(ncols)
    REAL(KIND=r8),    INTENT(IN   ) :: QCF(ncols,kMax)
    REAL(KIND=r8),    INTENT(IN   ) :: QCL(ncols,kMax)


    REAL(KIND=r8)    :: z0    (ncols) 

    REAL(KIND=r8)    :: gt    (ncols) 
    REAL(KIND=r8)    :: gq    (ncols) 
    REAL(KIND=r8)    :: speedm  (ncols)
    REAL(KIND=r8)    :: ah    (ncols)
    REAL(KIND=r8)    :: al    (ncols)
    REAL(KIND=r8)    :: am    (ncols)
    REAL(KIND=r8)    :: cuni  (ncols)
    REAL(KIND=r8)    :: cui   (ncols)
    REAL(KIND=r8)    :: cu    (ncols)
    REAL(KIND=r8)    :: ctni  (ncols)
    REAL(KIND=r8)    :: cti   (ncols)
    REAL(KIND=r8)    :: ct    (ncols)
    REAL(KIND=r8)    :: um    (ncols)
    REAL(KIND=r8)    :: vm    (ncols)
    REAL(KIND=r8)    :: tha   (ncols)
    REAL(KIND=r8)    :: thm   (ncols)
    REAL(KIND=r8)    :: dzm   (ncols)
    REAL(KIND=r8)    :: thvgm (ncols)
    REAL(KIND=r8)    :: rib   (ncols)
    REAL(KIND=r8)    :: ustar (ncols)
    REAL(KIND=r8)    :: gtsav (ncols)
    REAL(KIND=r8)    :: gqsav (ncols)
    REAL(KIND=r8)    :: tmsav (ncols)
    REAL(KIND=r8)    :: qmsav (ncols)
    REAL(KIND=r8)    :: tssav (ncols)
    REAL(KIND=r8)    :: dqg0  (ncols)
    REAL(KIND=r8)    :: b00   (ncols)
    REAL(KIND=r8)    :: b03   (ncols)
    REAL(KIND=r8)    :: b04   (ncols)
    REAL(KIND=r8)    :: c0    (ncols)
    REAL(KIND=r8)    :: b30   (ncols)
    REAL(KIND=r8)    :: b33   (ncols)
    REAL(KIND=r8)    :: c3    (ncols)
    REAL(KIND=r8)    :: b40   (ncols)
    REAL(KIND=r8)    :: b44   (ncols)
    REAL(KIND=r8)    :: c4    (ncols)
    REAL(KIND=r8)    :: rmi_uk  (ncols) 
    REAL(KIND=r8)    :: rhi_uk  (ncols) 
    REAL(KIND=r8)    :: evap_uk    (ncols) 
    REAL(KIND=r8)    :: sens_uk    (ncols)    
    REAL(KIND=r8)    :: ustar_uk(ncols) 
    REAL(KIND=r8)    :: RADNET(ncols)
    !LOGICAL :: jstneu
    !REAL(KIND=r8) :: u2 (ncols)

    INTEGER :: i
    INTEGER :: ncount
    REAL(KIND=r8)    :: gbyhl
    REAL(KIND=r8)    :: gbycp
    REAL(KIND=r8)    :: gr100
    REAL(KIND=r8)    :: gb100
    REAL(KIND=r8)    :: zgrav
    REAL(KIND=r8)    :: gice
    REAL(KIND=r8)    :: hscap
    REAL(KIND=r8)    :: sl1kap
    REAL(KIND=r8)    :: st4
    REAL(KIND=r8)    :: dti
    REAL(KIND=r8)    :: dtm
    REAL(KIND=r8)    :: dtmdt
    REAL(KIND=r8)    :: dqm
    REAL(KIND=r8)    :: dqmdt
    !*JPB REAL(KIND=r8), PARAMETER :: dd=0.05_r8
    REAL(KIND=r8), PARAMETER :: dd=3.0_r8 ! Total depth of the ice slab (m), Using ECMWF value
    REAL(KIND=r8), PARAMETER :: tice=271.16_r8
    REAL(KIND=r8), PARAMETER :: dice=2.0_r8
    REAL(KIND=r8), PARAMETER :: hice=2.03_r8
    REAL(KIND=r8), PARAMETER :: rhoice=920.0_r8 ! Mean ice density (kg/m3)
    REAL(KIND=r8), PARAMETER :: cice=2093.0_r8  ! Heat Capacity of Ice (J/Kg)


    sens =0.0_r8
    evap =0.0_r8
    rmi_uk  =0.0_r8
    rhi_uk  =0.0_r8
evap_uk=0.0_r8
sens_uk=0.0_r8
    ustar_uk=0.0_r8
    umom =0.0_r8
    vmom =0.0_r8
    rmi  =0.0_r8
    rhi  =0.0_r8
    rnet =0.0_r8
    z0   =0.0_r8
    qsfc =0.0_r8
    tsfc =0.0_r8
    bstar=0.0_r8
    gt    =0.0_r8
    gq    =0.0_r8 
    speedm=0.0_r8
    ah    =0.0_r8
    al    =0.0_r8
    am    =0.0_r8
    cuni  =0.0_r8
    cui   =0.0_r8
    cu    =0.0_r8
    ctni  =0.0_r8
    cti   =0.0_r8
    ct    =0.0_r8
    um    =0.0_r8
    vm    =0.0_r8
    tha   =0.0_r8
    thm   =0.0_r8
    dzm   =0.0_r8
    thvgm =0.0_r8
    rib   =0.0_r8
    ustar =0.0_r8
    gtsav =0.0_r8
    gqsav =0.0_r8
    tmsav =0.0_r8
    qmsav =0.0_r8
    tssav =0.0_r8
    dqg0  =0.0_r8
    b00   =0.0_r8
    b03   =0.0_r8
    b04   =0.0_r8
    c0    =0.0_r8
    b30   =0.0_r8
    b33   =0.0_r8
    c3    =0.0_r8
    b40   =0.0_r8
    b44   =0.0_r8
    c4    =0.0_r8
    Ustarm=0.0_r8
    rho   =0.0_r8

    gbyhl =0.0_r8
    gbycp =0.0_r8
    gr100 =0.0_r8
    gb100 =0.0_r8
    zgrav =0.0_r8
    gice =0.0_r8
    hscap =0.0_r8
    sl1kap =0.0_r8
    st4 =0.0_r8
    dti =0.0_r8
    dtm =0.0_r8
    dtmdt =0.0_r8
    dqm =0.0_r8
    dqmdt =0.0_r8

    gr100 =rgas*0.01_r8
    gbycp =grav/(cp*delsig*100.0_r8 *sigki)
    gbyhl =grav/(hltm*delsig*100.0_r8 )
    gb100 =grav/(   delsig*100.0_r8 )
    zgrav =0.032_r8 /grav
    gice  =hice/dice ! 2.03_r8/2.0_r8
    hscap =rhoice*cice*dd/dtc3x
    sl1kap=sigki
    st4   =stefan*4.0_r8
    dti   =1.0_r8 /dtc3x
    gt=t_sib(:,1)
    gq=sh_sib(:,1)
    DO i = 1, ncols
       RADNET(i)=ySwSfcNet(i)+LwSfcNet(i)
       IF(mskant(i) == 1_i8)THEN
          rnet (i)=-697.58_r8*slrad(i)
          rho  (i)=gps(i)/(gr100*gt(i))
          ah   (i)=gbycp/gps(i)
          al   (i)=gbyhl/gps(i)
          dqg0 (i)=0.622_r8 *EXP(30.25353_r8 -5418.0_r8 /tsurf(i)) &
               /(tsurf(i)*tsurf(i)*gps(i))
          gtsav(i)=gt   (i)
          gqsav(i)=gq   (i)
          tssav(i)=tsurf(i)
          tmsav(i)=tmtx (i,3)
          qmsav(i)=qmtx (i,3)
       END IF
    END DO

    c0  =0.0_r8
    cond=0.0_r8
    stor=0.0_r8

    ncount=0
8000 CONTINUE
    ncount=ncount+1
    !
    !     the first call to vntlat just gets the neutral values of ustar
    !     and ventmf.
    !

    !    jstneu=.TRUE.
    !    CALL vntlt2 &
    !       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
    !       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
    !       thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,jstneu,u2  )
    !
    !   jstneu=.FALSE.
    !    CALL vntlt2 &
    !       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
    !       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
    !       thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,jstneu,u2  )
    !
    CALL vntlt1_cola ( &
         rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
         sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
         thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant ,kMax)
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF(z0sea  (i)<= 0.0_r8) z0sea(i)=z0(i)
       END IF
    END DO


    CALL  SF_EXCH (&
                   nCols     , &
                   kMAx      , &
                   cldtot    , &
                   QCF       , &
                   QCL       , &
                   sh_sib    , &
                   T_sib     , &
                   gu        , &
                   gv        , &
                   sinclt    , &
                   gps       , &
                   radnet    , &
                   tsurf     , &
                   pblh      , &
                   mskant    , &
                   iMask     , &
                   ABS(tsurf), & !TSTAR_LAND , &
                   ABS(tsea) , & !TSTAR_SSI , &
                   ABS(tsea) , & !TSTAR_SEA ,&
                   ABS(tsea) , & !TSTAR_SICE ,&
                   rmi_uk    , & 
                   rhi_uk    , &
                   evap_uk , &
                   sens_uk , &
                   ustar_uk   , &
                   z0sea      & !Z0MSEA &
                               )

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt  (i)    =gtsav(i)
          gq  (i)    =gqsav(i)
          tsurf(i)   =tssav(i)
          tmtx(i,3)=tmsav(i)
          qmtx(i,3)=qmsav(i)
       END IF
    END DO
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             !
             ! Solution of sea ice
             !
             !
             !$$$$$$->PK
             rhi  (i) = rhi_uk(i)
             rmi  (i) = rmi_uk(i)
             !$$$$$$->PK
             !
             b00(i)=   hscap+cp*rho(i)*rhi(i) &
                  +hltm*rho(i)*rhi(i)*dqg0(i) &
                  +gice+st4*tsurf(i)**3
             b03(i)=        -cp*rho(i)*rhi(i)*sl1kap
             b04(i)=-hltm*rho(i)*rhi(i)
             !
             ! Right side of eq.41 section III.A 
             ! COLA Physics Description Manual
             !
             c0 (i)=rnet(i) -cp*rho(i)*rhi(i)*(tsurf(i)-sl1kap*gt(i)) &
                  -hltm*rho(i)*rhi(i)*(qsurf(i)-       gq(i)) &
                  -gice*(tsurf(i)-tice)-stefan*tsurf(i)**4
             b30(i)=               -ah (i)*cp*rho(i)*rhi(i)
             b33(i)=tmtx(i,2)*dti-b30(i)*          sl1kap
             c3 (i)=tmtx(i,3)    -b30(i)*(tsurf(i)-sl1kap*gt(i))
             b40(i)=               -al(i)*hltm*rho(i)*rhi(i)* dqg0 (i)
             b44(i)=qmtx(i,2)*dti+al(i)*hltm*rho(i)*rhi(i)
             c4 (i)=qmtx(i,3)    + &
                  al(i)*hltm*rho(i)*rhi(i)*(qsurf(i)-gq(i))
             b00(i)=b00(i)-b30(i)*b03(i)/b33(i)-b40(i)*b04(i)/b44(i)
             c0 (i)=c0 (i)-c3 (i)*b03(i)/b33(i)-c4 (i)*b04(i)/b44(i)
             c0 (i)=c0 (i)/b00(i)
             tsurf(i)=tsurf(i)+c0(i)
             tmtx(i,3)=(c3(i)-b30(i)*c0(i))/(b33(i)*dtc3x)
             qmtx(i,3)=(c4(i)-b40(i)*c0(i))/(b44(i)*dtc3x)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
             !
             ! Solution of sea water
             !
             !$$$$$$->PK
             rhi  (i) = rhi_uk(i)
             rmi  (i) = rmi_uk(i)
             !$$$$$$->PK
             zorl (i) = 100.0_r8    * zgrav*speedm(i)*rhi(i)
             !$$$$$$->PK
             sens (i) = sens_uk (i)
             evap (i) = evap_uk (i)    
             !pk  sens (i)= rho(i)*cp*(tsurf(i)-gt(i)*sigki(1))*rhi(i)
             !pk  evap (i)= rho(i)*hl*(qsurf(i)-gq(i))*rhi(i)
             !$$$$$$->PK
             !zorl (i)= 100.0_r8 *zgrav*speedm(i)*rhi(i)
             !sens (i)= rho(i)*cp*(tsurf(i)-gt(i)*sigki)*rhi(i)
             !evap (i)= rho(i)*hltm*(qsurf(i)-gq(i))*rhi(i)
             tmtx(i,3)=(tmtx(i,3)+ah(i)*sens(i)) &
                  /(tmtx(i,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i))
             qmtx(i,3)=(qmtx(i,3)+al(i)*evap(i)) &
                  /(qmtx(i,2)+dtc3x*al(i)*rho(i)*hltm*rhi(i))
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          gt(i)=gt(i)+tmtx(i,3)*dtc3x
          gq(i)=gq(i)+qmtx(i,3)*dtc3x
       END IF
    END DO

    IF (ncount == 1) go to 8000

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN          
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             !
             ! Solution of sea ice
             !
             rhi  (i) = rhi_uk(i)
             rmi  (i) = rmi_uk(i)
             sens (i) = sens_uk(i)
             evap (i) = evap_uk(i)

             !sens(i)=rho(i)*cp*  (tsurf(i)   -  gt(i)*sigki   )*rhi(i)
             !evap(i)=rho(i)*hltm*(qsurf(i)   -  gq(i)         )*rhi(i)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) > tice+0.01_r8) THEN
             !
             ! Solution of sea water
             !
             rhi  (i) = rhi_uk(i)
             rmi  (i) = rmi_uk(i)
             sens (i) = sens_uk(i)
             evap (i) = evap_uk(i)
          END IF       
          !sens(i)=rho(i)*cp*(tsurf(i)-gt(i)*sigki)*rhi(i)
          !evap(i)=rho(i)*hltm*(qsurf(i)-gq(i)         )*rhi(i)
          dtmdt=(ah(i)*sens(i))/(tmtx(i,2)+dtc3x*ah(i)*rho(i)*cp*rhi(i)+0.0001)
          dqmdt=(al(i)*evap(i))/(qmtx(i,2)+dtc3x*al(i)*rho(i)*hltm*rhi(i)+0.0001)
          dtm=dtmdt*dtc3x
          dqm=dqmdt*dtc3x
          tsfc   (i)=gt(i)+dtm
          qsfc   (i)=gq(i)+dqm

          gt  (i)=gtsav(i)
          gq  (i)=gqsav(i)
          rnet(i)=rnet(i)-stefan*tsurf(i)**4
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice+0.01_r8) THEN
             cond(i)=gice*(tsurf(i)-tice)
             stor(i)=hscap*c0(i)
             rnet(i)=rnet(i)-stefan*tsurf(i)**3*4.0_r8 *c0(i)
             tsurf(i)=MIN(tsurf(i),tice)
             tsea (i)=-   tsurf(i)
          END IF
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          !$$$$$$->PK
          !rmi  (i)= WSPD_SEA(i)*(CM_SEA(i) )
          !$$$$$$->PK
          umom(i)=rho(i)*gu(i,1)*rmi(i)
          vmom(i)=rho(i)*gv(i,1)*rmi(i)
          am  (i)=gb100/gps(i)
          umtx(i,3)=(umtx(i,3)-am(i)*umom(i)) &
               /(umtx(i,2)+dtc3x*am(i)*rho(i)*rmi(i))
          umtx(i,4)=(umtx(i,4)-am(i)*vmom(i)) &
               /(umtx(i,2)+dtc3x*am(i)*rho(i)*rmi(i))
          !
          !     set surface stress use of pseudo winds to true winds
          !     for output diagnostics
          !
          umom(i)=umom(i)/sinclt(i)
          vmom(i)=vmom(i)/sinclt(i)
          Ustarm(i) = SQRT(umom(i)**2 + vmom(i)**2)
          IF(Ustarm(i)==0.0_r8)Ustarm(i)=0.007_r8
          um  (i)=gu (i,1)/sinclt(i)
          vm  (i)=gv (i,1)/sinclt(i)
          speedm(i)=SQRT(um(i)**2 + vm(i)**2)
          speedm(i)=MAX(2.0_r8 , speedm(i))
          dzm   (i)=rbyg*gt(i)
          bstar(i)=cu(i)*grav*(ct(i)*(tsfc(i)-gt(i)-(grav/cp)*dzm(i))/gt(i))! + mapl_vireps*ct(i)*(qsfc(i)-gq(i)))


       END IF
    END DO
  END SUBROUTINE seasfc_ukme


  ! vntlt1 :performs ventilation mass flux, based on deardorff, mwr, 1972?.



  SUBROUTINE vntlt1_wgfs &
       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
       thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,kMax  )
    !
    !==========================================================================
    !==========================================================================
    !==========================================================================
    ! imax..........number of grid points on a gaussian latitude circle
    ! z0ice.........Roughness length of ice
    ! sinclt........sinclt=SIN(colrad(latitu))
    ! rmi...........rmi   (i)=cu(i)*ustar(i), where
    !               cu is friction  transfer coefficients
    !               ustar is surface friction velocity  (m/s)
    ! rhi...........rhi   (i)=ct(i)*ustar(i), where
    !               ct is heat transfer coefficients.
    !               ustar is surface friction velocity  (m/s)
    ! gu............(zonal      velocity)*sin(colat)
    ! gv............(meridional velocity)*sin(colat)
    ! gt............temperature
    ! tsurf.........earth's surface temperature used for radiation
    !               for the first time step when ground temperature is not yet
    !               computed (this is done by subr.tsinit ),
    ! tsea..........effective surface radiative temperature ( tgeff )
    ! zorl..........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !               zgrav =0.032 /grav
    ! delsig
    ! sigki ........sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !               sigma coordinate at middle of layer and akappa=gasr/cp
    ! cuni..........neutral friction transfer  coefficients.
    ! cui...........cui   (i)=cuni(i)*EXP( aa-SQRT(aa*aa+bb*f))
    !               cui   (i)=cuni(i)*EXP(-tt+SQRT(tt*tt+ss*f))
    ! cu............Friction  transfer coefficients.
    ! ctni..........neutral heat transfer coefficients.
    ! cti...........cti   (i)=ctni(i)*EXP( qq-SQRT(qq*qq+rr*g))
    !               cti   (i)=cui (i)
    ! ct............heat transfer coefficients.
    ! speedm........speedm(i)=SQRT(gu(i)**2+gv(i)**2)*sincli, where
    !               sincli=1.0 /sinclt
    ! tha...........tha   (i)= tsurf(i)
    ! thm...........thm   (i)= gt(i)*sigki(1)
    ! dzm...........dzm   (i)=gt(i)*rbyg
    !               rbyg  =gasr/grav*delsig(1)*0.5
    ! thvgm.........thvgm (i)= tha(i)-thm(i)
    ! rib...........bulk richardson number.
    ! z0............Roughness length
    ! ustarr........surface friction velocity  (m/s)
    ! gasr..........gas constant of dry air        (j/kg/k)
    ! grav..........grav   gravity constant        (m/s**2)
    !==========================================================================
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rmi   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rhi   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gu    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gv    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gt    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsurf (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsea  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: zorl  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: sigki  
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cui   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cti   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: speedm(ncols)
    REAL(KIND=r8),    INTENT(inout) :: tha   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thvgm (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rib   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    INTEGER(KIND=i8), INTENT(in   ) :: mskant(ncols)

    REAL(KIND=r8),        PARAMETER :: vkrmn=0.40_r8
    REAL(KIND=r8),        PARAMETER :: ribc=3.05_r8
    REAL(KIND=r8),        PARAMETER :: aa=1.2270_r8
    REAL(KIND=r8),        PARAMETER :: bb=1.2642_r8
    REAL(KIND=r8),        PARAMETER :: tt=1.8900_r8
    REAL(KIND=r8),        PARAMETER :: ss=5.0519_r8
    REAL(KIND=r8),        PARAMETER :: ee=1.2743_r8
    REAL(KIND=r8),        PARAMETER :: ff=3.4805_r8
    REAL(KIND=r8),        PARAMETER :: gg=0.87581_r8
    REAL(KIND=r8),        PARAMETER :: hh=-1.5630_r8
    REAL(KIND=r8),        PARAMETER :: pp=10.815_r8
    REAL(KIND=r8),        PARAMETER :: qq=1.3462_r8
    REAL(KIND=r8),        PARAMETER :: rr=1.8380_r8
    REAL(KIND=r8)                   :: sincli(ncols)
    REAL(KIND=r8)                   :: f
    REAL(KIND=r8)                   :: g
    INTEGER                :: i

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          z0(i)=0.001_r8
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) >= 271.17_r8) THEN
             z0(i)=MAX(0.01_r8*zorl(i),1.0e-3_r8    )
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < 271.17_r8) THEN
             z0(i)=z0ice
          END IF
          sincli(i)=1.0_r8 /sinclt(i)
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) <= 0.0_r8) THEN
          speedm(i)=SQRT(gu(i,1)**2+gv(i,1)**2)*sincli(i)
          speedm(i)=MAX(2.0_r8 ,speedm(i))
          dzm   (i)=gt(i)*rbyg
          cuni(i)=LOG(dzm(i)/z0(i))/vkrmn*gg+hh
          ctni(i)=cuni(i)
          !
          !     stability branch based on bulk richardson number.
          !
          thm   (i)= gt(i)*sigki
          tha   (i)= tsurf(i)
          thvgm (i)= tha(i)-thm(i)
          rib   (i)=-thvgm(i)*grav*dzm(i)/ (thm(i)*speedm(i)**2)
          rib   (i)=MAX(-1.25_r8 ,rib(i))
          rib   (i)=MIN( 1.25_r8 ,rib(i))
          IF (rib(i) < 0.0_r8) THEN
             f        =LOG(1.0_r8-ee*rib(i))
             cui   (i)=cuni(i)*EXP( aa-SQRT(aa*aa+bb*f))
             g        =LOG(1.0_r8-ff*rib(i))
             cti   (i)=ctni(i)*EXP( qq-SQRT(qq*qq+rr*g))
          ELSE
             f        =LOG(1.0_r8+pp*rib(i))
             cui   (i)=cuni(i)*EXP(-tt+SQRT(tt*tt+ss*f))
             cti   (i)=cui (i)
          END IF
          cu    (i)=1.0_r8/cui(i)
          ct    (i)=1.0_r8/cti(i)
          !
          !     surface friction velocity and ventilation mass flux
          !
          ustar (i)=speedm(i)*cu(i)
          rmi   (i)=cu(i)*ustar(i)
          rhi   (i)=ct(i)*ustar(i)
          END IF
       END IF
    END DO
  END SUBROUTINE vntlt1_wgfs



  ! vntlt1 :performs ventilation mass flux, based on deardorff, mwr, 1972?.



  SUBROUTINE vntlt1_cola &
       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
       thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,kMax  )
    !
    !==========================================================================
    !==========================================================================
    !==========================================================================
    ! imax..........number of grid points on a gaussian latitude circle
    ! z0ice.........Roughness length of ice
    ! sinclt........sinclt=SIN(colrad(latitu))
    ! rmi...........rmi   (i)=cu(i)*ustar(i), where
    !               cu is friction  transfer coefficients
    !               ustar is surface friction velocity  (m/s)
    ! rhi...........rhi   (i)=ct(i)*ustar(i), where
    !               ct is heat transfer coefficients.
    !               ustar is surface friction velocity  (m/s)
    ! gu............(zonal      velocity)*sin(colat)
    ! gv............(meridional velocity)*sin(colat)
    ! gt............temperature
    ! tsurf.........earth's surface temperature used for radiation
    !               for the first time step when ground temperature is not yet
    !               computed (this is done by subr.tsinit ),
    ! tsea..........effective surface radiative temperature ( tgeff )
    ! zorl..........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !               zgrav =0.032 /grav
    ! delsig
    ! sigki ........sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !               sigma coordinate at middle of layer and akappa=gasr/cp
    ! cuni..........neutral friction transfer  coefficients.
    ! cui...........cui   (i)=cuni(i)*EXP( aa-SQRT(aa*aa+bb*f))
    !               cui   (i)=cuni(i)*EXP(-tt+SQRT(tt*tt+ss*f))
    ! cu............Friction  transfer coefficients.
    ! ctni..........neutral heat transfer coefficients.
    ! cti...........cti   (i)=ctni(i)*EXP( qq-SQRT(qq*qq+rr*g))
    !               cti   (i)=cui (i)
    ! ct............heat transfer coefficients.
    ! speedm........speedm(i)=SQRT(gu(i)**2+gv(i)**2)*sincli, where
    !               sincli=1.0 /sinclt
    ! tha...........tha   (i)= tsurf(i)
    ! thm...........thm   (i)= gt(i)*sigki(1)
    ! dzm...........dzm   (i)=gt(i)*rbyg
    !               rbyg  =gasr/grav*delsig(1)*0.5
    ! thvgm.........thvgm (i)= tha(i)-thm(i)
    ! rib...........bulk richardson number.
    ! z0............Roughness length
    ! ustarr........surface friction velocity  (m/s)
    ! gasr..........gas constant of dry air        (j/kg/k)
    ! grav..........grav   gravity constant        (m/s**2)
    !==========================================================================
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kMax
    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rmi   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rhi   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gu    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gv    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gt    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsurf (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsea  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: zorl  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: sigki  
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cui   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cti   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: speedm(ncols)
    REAL(KIND=r8),    INTENT(inout) :: tha   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thvgm (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rib   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    INTEGER(KIND=i8), INTENT(in   ) :: mskant(ncols)

    REAL(KIND=r8),        PARAMETER :: vkrmn=0.40_r8
    REAL(KIND=r8),        PARAMETER :: ribc=3.05_r8
    REAL(KIND=r8),        PARAMETER :: aa=1.2270_r8
    REAL(KIND=r8),        PARAMETER :: bb=1.2642_r8
    REAL(KIND=r8),        PARAMETER :: tt=1.8900_r8
    REAL(KIND=r8),        PARAMETER :: ss=5.0519_r8
    REAL(KIND=r8),        PARAMETER :: ee=1.2743_r8
    REAL(KIND=r8),        PARAMETER :: ff=3.4805_r8
    REAL(KIND=r8),        PARAMETER :: gg=0.87581_r8
    REAL(KIND=r8),        PARAMETER :: hh=-1.5630_r8
    REAL(KIND=r8),        PARAMETER :: pp=10.815_r8
    REAL(KIND=r8),        PARAMETER :: qq=1.3462_r8
    REAL(KIND=r8),        PARAMETER :: rr=1.8380_r8
    REAL(KIND=r8)                   :: sincli(ncols)
    REAL(KIND=r8)                   :: f
    REAL(KIND=r8)                   :: g
    INTEGER                :: i

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          z0(i)=0.001_r8
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) >= 271.17_r8) THEN
             z0(i)=0.01_r8*zorl(i)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < 271.17_r8) THEN
             z0(i)=z0ice
          END IF
          sincli(i)=1.0_r8 /sinclt(i)
       END IF
    END DO

    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) <= 0.0_r8) THEN
          speedm(i)=SQRT(gu(i,1)**2+gv(i,1)**2)*sincli(i)
          speedm(i)=MAX(2.0_r8 ,speedm(i))
          dzm   (i)=gt(i)*rbyg
          cuni(i)=LOG(dzm(i)/z0(i))/vkrmn*gg+hh
          ctni(i)=cuni(i)
          !
          !     stability branch based on bulk richardson number.
          !
          thm   (i)= gt(i)*sigki
          tha   (i)= tsurf(i)
          thvgm (i)= tha(i)-thm(i)
          rib   (i)=-thvgm(i)*grav*dzm(i)/ (thm(i)*speedm(i)**2)
          rib   (i)=MAX(-1.25_r8 ,rib(i))
          rib   (i)=MIN( 1.25_r8 ,rib(i))
          IF (rib(i) < 0.0_r8) THEN
             f        =LOG(1.0_r8-ee*rib(i))
             cui   (i)=cuni(i)*EXP( aa-SQRT(aa*aa+bb*f))
             g        =LOG(1.0_r8-ff*rib(i))
             cti   (i)=ctni(i)*EXP( qq-SQRT(qq*qq+rr*g))
          ELSE
             f        =LOG(1.0_r8+pp*rib(i))
             cui   (i)=cuni(i)*EXP(-tt+SQRT(tt*tt+ss*f))
             cti   (i)=cui (i)
          END IF
          cu    (i)=1.0_r8/cui(i)
          ct    (i)=1.0_r8/cti(i)
          !
          !     surface friction velocity and ventilation mass flux
          !
          ustar (i)=speedm(i)*cu(i)
          rmi   (i)=cu(i)*ustar(i)
          rhi   (i)=ct(i)*ustar(i)
          END IF
       END IF
    END DO
  END SUBROUTINE vntlt1_cola





  ! vntlt1 :performs ventilation mass flux, based on deardorff, mwr, 1972?.



  SUBROUTINE vntlt2 &
       (rmi   ,rhi   ,gu    ,gv    ,gt    ,tsurf ,tsea  ,ncols , &
       sigki ,cuni  ,cui   ,cu    ,ctni  ,cti   ,ct    ,speedm,tha   , &
       thm   ,dzm   ,thvgm ,rib   ,z0    ,zorl  ,ustar ,sinclt,mskant,jstneu,u2 ,kMax )
    !
    !==========================================================================
    !==========================================================================
    !==========================================================================
    ! imax..........number of grid points on a gaussian latitude circle
    ! z0ice.........Roughness length of ice
    ! sinclt........sinclt=SIN(colrad(latitu))
    ! rmi...........rmi   (i)=cu(i)*ustar(i), where
    !               cu is friction  transfer coefficients
    !               ustar is surface friction velocity  (m/s)
    ! rhi...........rhi   (i)=ct(i)*ustar(i), where
    !               ct is heat transfer coefficients.
    !               ustar is surface friction velocity  (m/s)
    ! gu............(zonal      velocity)*sin(colat)
    ! gv............(meridional velocity)*sin(colat)
    ! gt............temperature
    ! tsurf.........earth's surface temperature used for radiation
    !               for the first time step when ground temperature is not yet
    !               computed (this is done by subr.tsinit ),
    ! tsea..........effective surface radiative temperature ( tgeff )
    ! zorl..........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
    !               zgrav =0.032 /grav
    ! delsig
    ! sigki ........sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
    !               sigma coordinate at middle of layer and akappa=gasr/cp
    ! cuni..........neutral friction transfer  coefficients.
    ! cui...........cui   (i)=cuni(i)*EXP( aa-SQRT(aa*aa+bb*f))
    !               cui   (i)=cuni(i)*EXP(-tt+SQRT(tt*tt+ss*f))
    ! cu............Friction  transfer coefficients.
    ! ctni..........neutral heat transfer coefficients.
    ! cti...........cti   (i)=ctni(i)*EXP( qq-SQRT(qq*qq+rr*g))
    !               cti   (i)=cui (i)
    ! ct............heat transfer coefficients.
    ! speedm........speedm(i)=SQRT(gu(i)**2+gv(i)**2)*sincli, where
    !               sincli=1.0 /sinclt
    ! tha...........tha   (i)= tsurf(i)
    ! thm...........thm   (i)= gt(i)*sigki(1)
    ! dzm...........dzm   (i)=gt(i)*rbyg
    !               rbyg  =gasr/grav*delsig(1)*0.5
    ! thvgm.........thvgm (i)= tha(i)-thm(i)
    ! rib...........bulk richardson number.
    ! z0............Roughness length
    ! ustarr........surface friction velocity  (m/s)
    ! gasr..........gas constant of dry air        (j/kg/k)
    ! grav..........grav   gravity constant        (m/s**2)
    !==========================================================================
    !
    INTEGER, INTENT(in   ) :: ncols
    INTEGER, INTENT(in   ) :: kMax
    INTEGER, PARAMETER :: r8 = 8
    INTEGER, PARAMETER :: i8 = 8

    REAL(KIND=r8),    INTENT(in   ) :: sinclt(ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rmi   (ncols)
    REAL(KIND=r8),    INTENT(inout  ) :: rhi   (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: gu    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gv    (ncols,kMax)
    REAL(KIND=r8),    INTENT(in   ) :: gt    (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsurf (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: tsea  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: zorl  (ncols)
    REAL(KIND=r8),    INTENT(in   ) :: sigki  
    REAL(KIND=r8),    INTENT(inout) :: cuni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cui   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cu    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ctni  (ncols)
    REAL(KIND=r8),    INTENT(inout) :: cti   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ct    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: speedm(ncols)
    REAL(KIND=r8),    INTENT(inout) :: tha   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: dzm   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: thvgm (ncols)
    REAL(KIND=r8),    INTENT(inout) :: rib   (ncols)
    REAL(KIND=r8),    INTENT(inout) :: z0    (ncols)
    REAL(KIND=r8),    INTENT(inout) :: ustar (ncols)
    INTEGER(KIND=i8)      , INTENT(in   ) :: mskant(ncols)
    LOGICAL, INTENT(in   ) :: jstneu
    REAL(KIND=r8),    INTENT(inout) :: u2 (ncols)

    REAL(KIND=r8),        PARAMETER :: vkrmn=0.40_r8
    REAL(KIND=r8),        PARAMETER :: ribc=3.05_r8
    REAL(KIND=r8),        PARAMETER :: aa=1.2270_r8
    REAL(KIND=r8),        PARAMETER :: bb=1.2642_r8
    REAL(KIND=r8),        PARAMETER :: tt=1.8900_r8
    REAL(KIND=r8),        PARAMETER :: ss=5.0519_r8
    REAL(KIND=r8),        PARAMETER :: ee=1.2743_r8
    REAL(KIND=r8),        PARAMETER :: ff=3.4805_r8
    REAL(KIND=r8),        PARAMETER :: gg=0.87581_r8
    REAL(KIND=r8),        PARAMETER :: hh=-1.5630_r8
    REAL(KIND=r8),        PARAMETER :: pp=10.815_r8
    REAL(KIND=r8),        PARAMETER :: qq=1.3462_r8
    REAL(KIND=r8),        PARAMETER :: rr=1.8380_r8
    REAL(KIND=r8), PARAMETER ::  fsc=66.85_r8
    REAL(KIND=r8), PARAMETER ::  ftc=0.904_r8
    REAL(KIND=r8), PARAMETER ::  fvc=0.315_r8

    REAL(KIND=r8)                   :: sincli(ncols)
    !REAL(KIND=r8)                   :: f
    !REAL(KIND=r8)                   :: g
    REAL(KIND=r8)                   :: zl
    REAL(KIND=r8)                   :: xct1
    REAL(KIND=r8)                   :: xct2
    REAL(KIND=r8)                   :: xctu1
    REAL(KIND=r8)                   :: xctu2
    REAL(KIND=r8)                   :: grib
    REAL(KIND=r8)                   :: grzl
    REAL(KIND=r8)                   :: grz2
    REAL(KIND=r8)                   :: fvv
    REAL(KIND=r8)                   :: ftt
    REAL(KIND=r8)                   :: rzl
    REAL(KIND=r8)                   :: rz2
    INTEGER                :: i
    REAL(KIND=r8) :: rfac
    REAL(KIND=r8) :: vkrmni
    REAL(KIND=r8) :: g2
    REAL(KIND=r8) :: z2(ncols)
    REAL(KIND=r8) :: d(ncols)
    REAL(KIND=r8) ::ustarn(ncols)
    rfac  =1.0e2_r8 /rair
    vkrmni=1.0_r8  /vkrmn
    g2 = 0.75_r8
    DO i = 1, ncols
       z2(i)  = 0.500_r8
       d (i)  = (0.500_r8+0.100_r8)/2.0_r8
       IF(mskant(i) == 1_i8)THEN
          z0(i)=0.001_r8
          IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) >= 271.17_r8) THEN
             z0(i)=0.01_r8*zorl(i)
          ELSE IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < 271.17_r8) THEN
             z0(i)=z0ice
          END IF
          sincli(i)=1.0_r8 /sinclt(i)
       END IF
    END DO



    IF (jstneu) THEN
       DO i = 1, ncols
          IF(mskant(i) == 1_i8)THEN
             speedm(i)=SQRT(gu(i,1)**2+gv(i,1)**2)*sincli(i)
             speedm(i)=MAX(0.5_r8 ,speedm(i))
             dzm   (i)=gt(i)*rbyg
             zl = z2(i) + 11.785_r8  * z0(i)
             cuni(i)=LOG((dzm(i)-d(i))/z0(i))*vkrmni
             ustarn(i)=speedm(i)/cuni(i)
             IF (zl < dzm(i)) THEN
                xct1 = LOG((dzm(i)-d(i))/(zl-d(i)))
                xct2 = LOG((zl-d(i))/z0(i))
                xctu1 = xct1
                xctu2 = LOG((zl-d(i))/(z2(i)-d(i)))
                ctni(i) = (xct1 + g2 * xct2) *vkrmni
             ELSE
                xct2 =  LOG((dzm(i)-d(i))/z0(i))
                xctu1 =  0.0_r8
                xctu2 =  LOG((dzm(i)-d(i))/(z2(i)-d(i)))
                ctni(i) = g2 * xct2 *vkrmni
             END IF
             !
             !     neutral values of ustar and ventmf
             !
             u2(i) = speedm(i) - ustarn(i)*vkrmni*(xctu1 + g2*xctu2)
          END IF
       END DO
       RETURN
    END IF
    !
    !     stability branch based on bulk richardson number.
    !
    DO i = 1, ncols
       IF(mskant(i) == 1_i8)THEN
          IF (tsea(i) <= 0.0_r8) THEN

             !
             !     freelm(i)=.false.
             !
             speedm(i)=SQRT(gu(i,1)**2+gv(i,1)**2)*sincli(i)
             speedm(i)=MAX(0.5_r8 ,speedm(i))
             dzm        (i)=gt(i)*rbyg
             thm        (i)= gt(i)*sigki
             tha        (i)= tsurf(i)
             thvgm (i)= tha(i)-thm(i)        
             zl       = z2(i) + 11.785_r8  * z0(i)
             rib(i)   =-thvgm(i)   *grav*(dzm(i)-d(i)) &
                  /( thm(i)*(speedm(i)-u2(i))**2)
             ! Manzi Suggestion:
             ! rib   (i)=max(-10.0_r8  ,rib(i))
             rib(i)      =MAX(-1.5_r8  ,rib(i)   )
             rib(i)      =MIN( 0.165_r8  ,rib (i)  )
             IF (rib(i)    < 0.0_r8) THEN
                grib = -rib(i)
                grzl = -rib(i)   * (zl-d(i))/(dzm(i)-d(i))
                grz2 = -rib(i)   * z0(i)/(dzm(i)-d(i))
                fvv = fvc*grib
                IF (zl < dzm(i)) THEN
                   ftt = (ftc*grib) + (g2-1.0_r8) * (ftc*grzl) - g2 * (ftc*grz2)
                ELSE
                   ftt = g2*((ftc*grib) - (ftc*grz2))
                END IF
                cui(i)   = cuni(i) - fvv
                cti(i)    = ctni(i) - ftt
             ELSE
                rzl = rib(i)   /(dzm(i)-d(i))*(zl-d(i))
                rz2 = rib(i)   /(dzm(i)-d(i))*z0(i)
                fvv = fsc*rib(i)
                IF (zl < dzm(i)) THEN
                   ftt = (fsc*rib(i)) + (g2-1) * (fsc*rzl) - g2 * (fsc*rz2)
                ELSE
                   ftt = g2 * ((fsc*rib(i)) - (fsc*rz2))
                END IF
                cui(i)    = cuni(i) + fvv
                cti(i)    = ctni(i) + ftt
             ENDIF
             cu    (i)=1.0_r8/cui(i)
             !**(JP)** ct is not used anywhere else
             ct    (i)=1.0_r8/cti(i)
             !
             !
             !     surface friction velocity and ventilation mass flux
             !
             ustar (i)=speedm(i)*cu(i)
             !**(JP)** ran is not used anywhere else
             !ran(i) = ctni(i) / ustarn(i)
             !ran(i) = MAX(ran(i), 0.8_r8 )
             rmi        (i)=cu(i)*ustar(i)
             rhi        (i)=ct(i)*ustar(i)
          END IF
       END IF
    END DO
  END SUBROUTINE vntlt2


  SUBROUTINE Albedo_IBIS(jb,nCols,zenith,tsea,imask,avisb       ,avisd     ,anirb      , &
       anird )

    IMPLICIT NONE
    INTEGER , INTENT(IN) :: jb
    INTEGER , INTENT(IN) :: nCols
    REAL(KIND=r8) , INTENT(IN   ) :: zenith (nCols)       ! cosine of solar zenith angle   
    REAL(KIND=r8) , INTENT(IN   ) :: tsea   (nCols)       ! cosine of solar zenith angle   
    INTEGER(KIND=i8),INTENT(IN   ) :: imask (ncols)
    REAL(KIND=r8),INTENT(OUT  ) :: avisb (ncols)
    REAL(KIND=r8),INTENT(OUT  ) :: avisd (ncols)
    REAL(KIND=r8),INTENT(OUT  ) :: anirb (ncols)
    REAL(KIND=r8),INTENT(OUT  ) :: anird (ncols)

    INTEGER       :: nsol         ! number of points in indsol
    REAL(KIND=r8) :: solu   (nCols)! solar flux (direct + diffuse) absorbed by upper canopy leaves per unit canopy area (W m-2)
    REAL(KIND=r8) :: sols   (nCols)! solar flux (direct + diffuse) absorbed by upper canopy stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: soll   (nCols)! solar flux (direct + diffuse) absorbed by lower canopy leaves and stems per unit canopy area (W m-2)
    REAL(KIND=r8) :: solg   (nCols)! solar flux (direct + diffuse) absorbed by unit snow-free soil (W m-2)
    REAL(KIND=r8) :: soli   (nCols)! solar flux (direct + diffuse) absorbed by unit snow surface (W m-2)
    REAL(KIND=r8) :: scalcoefl(nCols,4)   ! term needed in lower canopy scaling
    REAL(KIND=r8) :: scalcoefu(nCols,4)   ! term needed in upper canopy scaling
    INTEGER       :: indsol (nCols)         ! index of current strip for points with positive coszen
    REAL(KIND=r8) :: topparu(nCols)        ! total photosynthetically active raditaion absorbed by top leaves of upper canopy (W m-2)
    REAL(KIND=r8) :: topparl(nCols)        ! total photosynthetically active raditaion absorbed by top leaves of lower canopy (W m-2)
    REAL(KIND=r8) :: albsod (nCols)          ! direct  albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsoi (nCols)          ! diffuse albedo for soil surface (visible or IR)
    REAL(KIND=r8) :: albsnd (nCols)          ! direct  albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: albsni (nCols)          ! diffuse albedo for snow surface (visible or IR)
    REAL(KIND=r8) :: relod  (nCols)         ! upward direct radiation per unit icident direct beam on lower canopy (W m-2)
    REAL(KIND=r8) :: reloi  (nCols)         ! upward diffuse radiation per unit incident diffuse radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: reupd  (nCols)         ! upward direct radiation per unit incident direct radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: reupi  (nCols)         ! upward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: ablod  (nCols)          ! fraction of direct  radiation absorbed by lower canopy
    REAL(KIND=r8) :: abloi  (nCols)          ! fraction of diffuse radiation absorbed by lower canopy
    REAL(KIND=r8) :: flodd  (nCols)          ! downward direct radiation per unit incident direct radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: dummy  (nCols)          ! placeholder, always = 0: no direct flux produced for diffuse incident
    REAL(KIND=r8) :: flodi  (nCols)         ! downward diffuse radiation per unit incident direct radiation on lower canopy (W m-2)
    REAL(KIND=r8) :: floii  (nCols)         ! downward diffuse radiation per unit incident diffuse radiation on lower canopy
    REAL(KIND=r8) :: terml  (nCols,7)          ! term needed in lower canopy scaling
    REAL(KIND=r8) :: termu  (nCols,7)          ! term needed in upper canopy scaling
    REAL(KIND=r8) :: abupd  (nCols)       ! fraction of direct  radiation absorbed by upper canopy
    REAL(KIND=r8) :: abupi  (nCols)         ! fraction of diffuse radiation absorbed by upper canopy
    REAL(KIND=r8) :: fupdd  (nCols)         ! downward direct radiation per unit incident direct beam on upper canopy (W m-2)
    REAL(KIND=r8) :: fupdi  (nCols)         ! downward diffuse radiation per unit icident direct radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: fupii  (nCols)         ! downward diffuse radiation per unit incident diffuse radiation on upper canopy (W m-2)
    REAL(KIND=r8) :: fwetu  (nCols)         ! fraction of upper canopy leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqu  (nCols)         ! proportion of fwetu due to liquid
    REAL(KIND=r8) :: fwets  (nCols)         ! fraction of upper canopy stem area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliqs  (nCols)         ! proportion of fwets due to liquid
    REAL(KIND=r8) :: fwetl  (nCols)         ! fraction of lower canopy stem & leaf area wetted by intercepted liquid and/or snow
    REAL(KIND=r8) :: rliql  (nCols)         ! proportion of fwetl due to liquid
    REAL(KIND=r8) :: coszen (nCols)       ! cosine of solar zenith angle
    REAL(KIND=r8) :: f
    REAL(KIND=r8) :: ocealb
    INTEGER       :: ncount
    INTEGER       :: ib  ,i, npoi
    npoi=0
    DO i=1,nCols
       IF (iMask(i) >= 1_i8) THEN
          npoi=npoi+1 
          coszen(npoi) = zenith(i)
       END IF
    END DO

    !
    ! calculate areal fractions wetted by intercepted h2o
    !
    CALL fwetcal(npoi               , &! INTENT(IN   )
         fwetu (1:npoi)     , &! INTENT(OUT  )
         rliqu (1:npoi)     , &! INTENT(OUT  )
         fwets (1:npoi)     , &! INTENT(OUT  )
         rliqs (1:npoi)     , &! INTENT(OUT  )
         fwetl (1:npoi)     , &! INTENT(OUT  )
         rliql (1:npoi)     , &! INTENT(OUT  )
         wliqu (1:npoi,jb)  , &! INTENT(IN   )
         wliqumax           , &! INTENT(IN   ) ::
         wsnou (1:npoi,jb)  , &! INTENT(IN   ) ::
         wsnoumax           , &! INTENT(IN   ) ::
         tu   (1:npoi,jb)   , &! INTENT(IN   )
         wliqs(1:npoi,jb)   , &! INTENT(IN   )
         wliqsmax           , &! INTENT(IN   )
         wsnos(1:npoi,jb)   , &! INTENT(IN   )
         wsnosmax           , &! INTENT(IN   )
         ts   (1:npoi,jb)   , &! INTENT(IN   )
         wliql(1:npoi,jb)   , &! INTENT(IN   )
         wliqlmax           , &! INTENT(IN   )
         wsnol(1:npoi,jb)   , &! INTENT(IN   )
         wsnolmax           , &! INTENT(IN   )
         tl   (1:npoi,jb)   , &! INTENT(IN   )
         epsilon            , &! INTENT(IN   )
         tmelt              )  ! INTENT(IN   )
    !
    ! set up for solar calculations
    !
    CALL solset(npoi                 , &! INTENT(IN   )
         nsol                 , &! INTENT(OUT  )
         nband                , &! INTENT(IN   )
         solu     (1:npoi)    , &! INTENT(OUT  )
         sols     (1:npoi)    , &! INTENT(OUT  )
         soll     (1:npoi)    , &! INTENT(OUT  )
         solg     (1:npoi)    , &! INTENT(OUT  )
         soli     (1:npoi)    , &! INTENT(OUT  )
         scalcoefl(1:npoi,1:4), &! INTENT(OUT  )
         scalcoefu(1:npoi,1:4), &! INTENT(OUT  )
         indsol   (1:npoi)    , &! INTENT(OUT  )
         topparu  (1:npoi)    , &! INTENT(OUT  )
         topparl  (1:npoi)    , &! INTENT(OUT  )
         asurd    (1:npoi,1:nband,jb) , &! INTENT(OUT  )
         asuri    (1:npoi,1:nband,jb) , &! INTENT(OUT  )
         coszen   (1:npoi))      ! INTENT(IN   )  
    !
    ! solar calculations for each waveband
    !
    DO  ib = 1, nband
       !
       ! solsur sets surface albedos for soil and snow
       ! solalb performs the albedo calculations
       ! solarf uses the unit-incident-flux results from solalb
       ! to obtain absorbed fluxes sol[u,s,l,g,i] and 
       ! incident pars sunp[u,l]
       !
       CALL solsur (ib                               , &! INTENT(IN   )
            tmelt                            , &! INTENT(IN   )
            nsol                             , &! INTENT(IN   )
            albsod (1:npoi)                  , &! INTENT(OUt  )
            albsoi (1:npoi)                  , &! INTENT(OUt  )
            albsnd (1:npoi)                  , &! INTENT(OUt  )
            albsni (1:npoi)                  , &! INTENT(OUt  )
            indsol (1:npoi)                  , &! INTENT(IN   )
            wsoi   (1:npoi,1:nsoilay,jb)     , &! INTENT(IN   )
            wisoi  (1:npoi,1:nsoilay,jb)     , &! INTENT(IN   )
            albsav (1:npoi,jb)               , &! INTENT(IN   )
            albsan (1:npoi,jb)               , &! INTENT(IN   )
            tsno   (1:npoi,1:nsnolay,jb)     , &! INTENT(IN   )
            coszen (1:npoi)                  , &! INTENT(IN   )
            npoi                             , &! INTENT(IN   )
            nsoilay                          , &! INTENT(IN   )
            nsnolay                            )! INTENT(IN   )

       CALL solalb (ib               , &! INTENT(IN   )
            relod (1:npoi)   , &! INTENT(OUT  )
            reloi (1:npoi)   , &! INTENT(OUT  )
            indsol(1:npoi)   , &! INTENT(IN   )
            reupd (1:npoi)   , &! INTENT(OUT  )
            reupi (1:npoi)   , &! INTENT(OUT  )
            albsnd(1:npoi)   , &! INTENT(IN   )
            albsni(1:npoi)   , &! INTENT(IN   )
            albsod(1:npoi)   , &! INTENT(IN   )
            albsoi(1:npoi)   , &! INTENT(IN   )
            fl    (1:npoi,jb), &! INTENT(IN   )
            fu    (1:npoi,jb), &! INTENT(IN   )
            fi    (1:npoi,jb), &! INTENT(IN   )
            asurd (1:npoi,1:nband,jb) , &! INTENT(INOUT)! local
            asuri (1:npoi,1:nband,jb) , &! INTENT(INOUT)! local
            npoi             , &! INTENT(IN   )
            nband            , &! INTENT(IN   )
            nsol             , &! INTENT(IN   )
            ablod (1:npoi)   , &! INTENT(OUT  )
            abloi (1:npoi)   , &! INTENT(OUT  )
            flodd (1:npoi)   , &! INTENT(OUT  )
            dummy (1:npoi)   , &! INTENT(OUT  )
            flodi (1:npoi)   , &! INTENT(OUT  )
            floii (1:npoi)   , &! INTENT(OUT  )
            coszen(1:npoi)   , &! INTENT(IN   )
            terml (1:npoi,1:7)  , &! INTENT(OUT  )
            termu (1:npoi,1:7)  , &! INTENT(OUT  )
            lai   (1:npoi,1:2,jb), &! INTENT(IN   )
            sai   (1:npoi,1:2,jb), &! INTENT(IN   )
            abupd (1:npoi)      , &! INTENT(OUT  )
            abupi (1:npoi)      , &! INTENT(OUT  )
            fupdd (1:npoi)      , &! INTENT(OUT  )
            fupdi (1:npoi)      , &! INTENT(OUT  )
            fupii (1:npoi)      , &! INTENT(OUT  )
            fwetl (1:npoi)      , &! INTENT(IN   )
            rliql (1:npoi)      , &! INTENT(IN   )
            rliqu (1:npoi)      , &! INTENT(IN   )
            rliqs (1:npoi)      , &! INTENT(IN   )
            fwetu (1:npoi)      , &! INTENT(IN   )
            fwets (1:npoi)      , &! INTENT(IN   )
            rhoveg(1:nband,1:2) , &! INTENT(IN   )
            tauveg(1:nband,1:2) , &! INTENT(IN   )
            orieh (1:2)         , &! INTENT(IN   )
            oriev (1:2)         , &! INTENT(IN   )
            tl    (1:npoi,jb)   , &! INTENT(IN   )
            ts    (1:npoi,jb)   , &! INTENT(IN   )
            tu    (1:npoi,jb)   , &! INTENT(IN   )
            pi                  , &! INTENT(IN   )
            tmelt               , &! INTENT(IN   )
            epsilon                )! INTENT(IN   )


    END DO

    ncount=0
    DO i=1,ncols
       IF(imask(i).GE.1_i8) THEN
          ncount=ncount+1
          avisb(i)=asuri(ncount,1,jb)                   !asurd  (npoi,nband)   ! local  ! direct albedo of surface system
          avisd(i)=asurd(ncount,1,jb)                   !asuri  (npoi,nband)   ! local  ! diffuse albedo of surface system 
          anirb(i)=asuri(ncount,2,jb)
          anird(i)=asurd(ncount,2,jb)
       ELSE IF(ABS(tsea(i)).GE.271.16e0_r8 +0.01e0_r8) THEN
          f=MAX(zenith(i),0.0e0_r8 )
          ocealb=0.12347e0_r8 +f*(0.34667e0_r8+f*(-1.7485e0_r8 + &
               f*(2.04630e0_r8 -0.74839e0_r8 *f)))
          avisb(i)=ocealb
          avisd(i)=oceald
          anirb(i)=ocealb
          anird(i)=oceald
       ELSE
          avisb(i)=icealv
          avisd(i)=icealv
          anirb(i)=icealn
          anird(i)=icealn
       END IF
    END DO

  END SUBROUTINE Albedo_IBIS
  !
  !***************************************************************************
  !                      (imonth,iday,iyear)
  REAL(KIND=r8) FUNCTION julday (imonth,iday,iyear,tod)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) :: imonth
    INTEGER, INTENT(IN   ) :: iday
    INTEGER, INTENT(IN   ) :: iyear
    REAL(KIND=r8)   , INTENT(IN   ) :: tod
    !
    ! compute the julian day from a normal date
    !
    julday= iday  &
         + MIN(1,MAX(0,imonth-1))*31  &
         + MIN(1,MAX(0,imonth-2))*(28+(1-MIN(1,MOD(iyear,4))))  &
         + MIN(1,MAX(0,imonth-3))*31  &
         + MIN(1,MAX(0,imonth-4))*30  &
         + MIN(1,MAX(0,imonth-5))*31  &
         + MIN(1,MAX(0,imonth-6))*30  &
         + MIN(1,MAX(0,imonth-7))*31  &
         + MIN(1,MAX(0,imonth-8))*31  &
         + MIN(1,MAX(0,imonth-9))*30  &
         + MIN(1,MAX(0,imonth-10))*31  &
         + MIN(1,MAX(0,imonth-11))*30  &
         + MIN(1,MAX(0,imonth-12))*31  &
         + tod/86400.0

  END FUNCTION julday

END MODULE Sfc_Ibis_Interface
