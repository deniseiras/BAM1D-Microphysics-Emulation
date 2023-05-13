 SUBROUTINE physcs &
      (gt    ,gq    ,gu    ,gv    ,gps   ,tc0   ,tg0   ,td0   ,w0    ,capac0, &
       tcm   ,tgm   ,tdm   ,wm    ,capacm,ppli  ,ppci  ,gyu   ,gyv   ,gtd   , &
       gqd   ,htr   ,clr   ,imask ,rvisb ,rvisd ,rnirb ,rnird ,dlwbot,tsea  , &
       colrad,sig   ,sigml ,delsig,ifday ,tod   ,avisb ,avisd ,anirb ,anird , &
       dswtop,rs    ,ulwtop,gl0   ,zorl  ,gtt   ,gqq   ,sheleg,tseam ,omg   , &
       htrc  ,rvisbc,rvisdc,rnirbc,rnirdc,dlwclr,uswtpc,rsclr ,ultclr,convts, &
       convcs,convbs,convc ,convt ,convb ,ustr  ,vstr  ,latco ,yvisb ,yvisd , &
       ynirb ,ynird ,yswtop,yvisbc,yvisdc,ynirbc,ynirdc,cldsav,ssib  ,ncols , &
       kmax  ,sigki ,xvisb ,xvisd ,xnirb ,xnird ,xswtop,xvisbc,xvisdc,xnirbc, &
       xnirdc,lonrad,ps    ,var   ,sens  ,evap  ,cos2d ,intcosz)

  USE Radiation, ONLY: &
       spmrad,         &
       radtim,         &
       rqvirt


  USE Surface, ONLY: &
       fysiks,       &
       seasfc,       &
       sextrp,       &
       radalb,       &
       wheat ,       &
       xcover,       &
       x0x   ,       &
       xd    ,       &
       xdc   ,       &
       xbc   ,       &  
       zlt   ,       &  
       green ,       &  
       chil  ,       &
       rstpar,       &
       topt  ,       &
       tll   ,       &
       tu    ,       &
       defac ,       &
       ph1   ,       &
       ph2   

   USE Diagnostics, ONLY: &
       updia    ,         &
       lgaus    ,         &
       gaus     ,         &
       nucf     ,         &
       iavrq    ,         &
       itavl    ,         &
       nuavl    ,         &
       lvavl    ,         &
       dodia    ,         &
       inavl    ,         &
       ixavl    ,         &
       iclcd    ,         &
       incf     ,         &
       ixcf     ,         &
       kravl    ,         &
       jrcf     ,         &
       krcf     ,         &
       ngaus    ,         &
       mgaus    ,         &
       ndrq     ,         &
       ncdg     ,         &
       jxavl    ,         &
       jxcdg    ,         &
       numx     ,         &
       ndshf    ,         &
       ndlhf    ,         &
       ndustr   ,         &
       ndvstr   ,         &
       nddtvd   ,         &
       ndduvd   ,         &
       nddvvd   ,         &
       nddqvd   ,         &
       ndtsb    ,         &
       ndcld    ,         &
       nddtlw   ,         &
       nddtsw   ,         &
       ndlwdb   ,         &
       ndlwub   ,         &
       ndlwut   ,         &
       ndswdt   ,         &
       ndswdb   ,         &
       ndswea   ,         &
       ndswe    ,         &
       ndswut   ,         &
       ndswub   ,         &
       ndlwnb   ,         &
       ndintl   ,         &
       ndrnof   ,         &
       ndldbc   ,         &
       ndlutc   ,         &
       ndsdbc   ,         &
       ndseac   ,         &
       ndswec   ,         &
       ndsutc   ,         &
       ndsubc   ,         &
       ndlnbc   ,         &
       ndta     ,         &
       ndea     ,         &
       ndlhbs   ,         &
       ndctd    ,         &
       ndcin    ,         &
       ndcst    ,         &
       ndccn    ,         &
       ndcsn    ,         &
       ndclw    ,         &
       ndemi    ,         &
       ndsod    ,         &
       ndrm     ,         &
       ndrarh   ,         &
       ndrb     ,         &
       ndrd     ,         &
       ndrg     ,         &
       ndrs     ,         &
       ndtxgw   ,         &
       ndtygw   ,         &
       ndugwd   ,         &
       ndvgwd 
       
   USE Constants, ONLY :     &
       ndavl    ,         &
       pie      ,        &
       cp       ,        &
       hl       ,        &
       gasr     ,        &
       grav     ,        &
       solcon   ,        &
       stefan   ,        &
       tf       ,        &
       epsfac   ,        &
       tice     ,        &
       oceald   ,        &
       icealn   ,        &
       icealv

   USE PlanBoundLayer, ONLY: &
        ympbl0,              &
        ympbl1,              &
	InitGwdd,            &
	Gwdd
	
   USE GridHistory, ONLY:       &
        IsGridHistoryOn,        &
        StoreGridHistory,       &
        StoreMaskedGridHistory, &
        ngrarh,                 &
        ngrm,                   &
        ngctd,                  &
        nglwh,                  &
        ngswh,                  &
        ngulwt,                 &
        ngdlwb,                 &
        ngcsz,                  &
        ngdswt,                 &
        ngsvb,                  &
        ngsvd,                  &
        ngsnb,                  &
        ngsnd,                  &
        ngavb,                  &
        ngavd,                  &
        nganb,                  &
        ngand,                  &
        ngcnd,                  &
        ngstr,                  &
        ngcld,                  &
        ngdrg,                  &
        ngrd1,                  &
        ngrd2,                  &
        ngrb,                   &
        ngrd,                   &
        ngrc,                   &
        ngrg,                   &
        ngrs,                   &
        ngea,                   &
        ngta,                   &
        nghc,                   &
        nghg,                   &
        ngect,                  &
        ngegt,                  &
        ngeci,                  &
        ngegi,                  &
        nglhbs,                 &
        ngchf,                  &
        ngghf,                  &
        ngrof,                  &
        ngdtvd,                 &
        ngduvd,                 &
        ngdvvd,                 &
        ngdqvd,                 &
        nguswt,                 &
        ngulwb,                 &
        ngustr,                 &
        ngvstr,                 &
        ngshf,                  &
        nglhf,                  &
        dogrh

   USE InputOutput, ONLY: &
       nfprt , &
       nferr , &
       ifprt 

   USE Options, ONLY : &
       varcut,&
       dogwd ,&
       mxrdcc,&
       lcnvl ,&
       lthncl,&
       cdhl  ,&
       istrt ,&
       intgr ,&
       first ,&
       co2val,&
       delt  ,&
       filta ,&
       nfin0 ,&
       nfin1 ,&
       initlz,&
       nfcnv0,&
       nfcldr,&
       irad  ,&
       iccon ,&
       swint ,&
       trint ,&
       yrl   ,&
       monl  ,&
       dtc3x ,&
       epsflt,&
       intg  ,&
       nls   ,&
       nlcs  ,&
       maxtid,&       
       dt    ,&
       idate ,&
       idatec,&
       kt    ,&
       ktm   ,&
       jdt   
       
  
  !
  ! $Author: cptec $
  ! $Date: 2000/01/10 22:27:56 $
  ! $Revision: 1.1 $
  !
  ! physcs :main subroutine for turbulence closure
  !         hashvadahn radiation coupled 3-d model
  !         p.sellers  sib
  !         gps is in mb
  !==========================================================================
  ! ncols.....Number of grid points on a gaussian latitude circle  
  ! jmax......Number of gaussian latitudes  
  ! kmax......Number of sigma levels  
  ! nls..... .Number of layers in the stratosphere.   
  ! nlcs......nlcs =   30  
  ! maxtid..../include/T062L28/restim.inc:      
  !           constant integer, parameter maxtid=131760  
  ! mgaus.....mgaus=	1
  ! ngaus.....ngaus=   21
  !            constant integer ,parameter npmx=315
  ! imx.......=ncols+1 or ncols+2   :this dimension instead of ncols
  !              is used in order to avoid bank conflict of memory
  !              access in fft computation and make it efficient. the
  !              choice of 1 or 2 depends on the number of banks and
  !              the declared type of grid variable (real*4,real*8)
  !              to be fourier transformed.
  !              cyber machine has the symptom.
  !              cray machine has no bank conflict, but the argument
  !              'imx' in subr. fft991 cannot be replaced by ncols      
  ! ityp......Numero das classes de solo 13 
  ! imon......Max. number of month at year (12)
  ! icg.......Parameter of the vegetation  (icg=1 top e icg=2 bottom )
  ! iwv.......Compriment de onda iwv=1=visivel, iwv=2=infravermelho
  !           proximo, iwv=3 infravermelho termal
  ! idp.......Parameter to the layers of soils idp=1->3
  ! ibd.......Condiction of vegetation ibd=1 green / ibd=2 
  ! ndavl.....Constant ndavl=   100 
  ! ndrq......Constant ndrq =   100 
  ! ncdg......Constant ncdg =    20
  ! jxavl.....Constant jxavl=   150 
  ! jxcdg.....Constant jxcdg=   40 
  ! numx......Constant numx = nunits-1, where nunits=230
  ! gt........Temperature    
  ! gq........Specific humidity    
  ! gu........(zonal      velocity)*sin(colat)    
  ! gv........(meridional velocity)*sin(colat)    
  ! gps.......Surface pressure in mb   
  ! tc0.......Temperatura da copa "dossel"(K)   modificada   
  ! tg0.......Temperatura da superficie do solo  (K)   modificada   
  ! td0.......Temperatura do solo profundo (K)   modificada   
  ! w0(id)....Grau de saturacao de umidade do solo id=1 na camada superficial
  ! w0(id)....Grau de saturacao de umidade do solo id=2 na camada de raizes
  ! w0(id)....Grau de saturacao de umidade do solo id=3 na camada de drenagem  
  ! capac0(iv).Agua interceptada iv=1 no dossel "water store capacity
  !             of leaves"(m)  modificada
  ! capac0(iv).Agua interceptada iv=2 na cobertura do solo (m)   modificada
  ! tcm........Temperatura da copa "dossel"(K)  
  ! tgm........Temperatura da superficie do solo  (K)   
  ! tdm........Temperatura do solo profundo (K)     
  ! wm    
  ! capacm.....Agua interceptada iv=2 na cobertura do solo (m)   
  ! ppli.......Precipitation rate ( large scale )       (mm/s)  
  ! ppci.......Precipitation rate ( cumulus )           (mm/s)  
  ! gyu........-(dv/dt) negative of tendency of v*cos(latitude)   
  ! gyv........(du/dt) tendency of zonal wind * cos(latitude)   
  ! gtd   
  ! gqd........Specific humidity   
  ! htr........Heating rate due to shrt wave radiation in deg/sec   
  ! clr........Cooling rate due to long wave radiation in deg/sec   
  ! imask......mascara continetal
  ! rvisb......visible beam cloudy skies (refer to downward surface fluxes) 
  ! rvisd......visible diffuse cloudy skies (refer to downward surface fluxes) 
  ! rnirb......near-ir beam cloudy skies (refer to downward surface fluxes) 
  ! rnird......near-ir diffuse cloudy skies (refer to downward surface fluxes) 
  ! dlwbot.....downward longwave radiation at the bottom in w/m**2
  ! tsea.......effective surface radiative temperature ( tgeff )  
  ! colrad.....colatitude  colrad=0-3.14 from np to sp in radians
  ! sig........sigma coordinate at middle of layer   
  ! sigml......sigma coordinate at bottom of layer 
  ! delsig      k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(2)
  !             k=3/2----sigml,ric,rf,km,kh,b,l -----------
  !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(1)
  !             k=1/2----sigml ----------------------------
  !
  ! istrt.......istrt = jdt =time step in getdia
  ! ifday.......model forecast day
  ! tod.........model forecast time of day in seconds   
  ! avisb.......visible beam surface albedo 
  ! avisd.......visible diffuse surface albedo 
  ! anirb.......near-ir beam surface albedo 
  ! anird.......near-ir diffuse surface albedo 
  ! uswtop......shortwave upward at top 
  ! dswtop......swinc....solar input at top of atmosphere
  ! rs..........net surface ir radiation in w/m**2     
  ! ulwtop......long wave flux at top of atmosphere in w/m**2
  ! alon........define constant alon=0.0 at subroutine gfidi.f90: 
  ! dt........time interval,usually =delt,but changes
  !             in nlnmi (dt=1.) and at dead start(delt/4,delt/2) 
  ! intgr.......intgr=2  time integration of surface physical variable 
  !                      is done by leap-frog implicit scheme. this 
  !                      conseves enegy and h2o.
  !             intgr=1  time integration of surface physical variable
  !                      is done by backward implicit scheme.
  ! gl0.........maximum mixing length l0 in blackerdar's formula
  !             l=k0*z/(1+k0*z/l0)  
  ! zorl........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
  !             zgrav =0.032 /grav and rhi   (i)=ct(i)*ustar(i), where
  !             ct is heat transfer coefficients.
  !             ustar is surface friction velocity  (m/s)
  !             speedm(i)=SQRT(gu(i)**2+gv(i)**2)*sincli, where
  !             sincli=1.0 /sinclt
  !
  ! gtt.........gtt =  gtmp(imx,kmax) input  : temperature (fourier).
  !                                   output : "s" as given by eq. (19)
  !                                            in noaa tech report nws 30  
  ! gqq.........gqq = gq(imx,kmax)     input : specific humidity (fourier).
  !                                   output : tendency of specific humidity
  !                                            without sources and sinks
  !                                            (fourier).   
  ! sheleg......snow amount in mm (equivalent water depth)
  ! tseam.......tseam = gtsea (ncols,jmax)  input(gsstcd) lt.0=sea surface temp.
  !                                                      gt.0=ground temp. 
  ! omg.........omg   =  vertical velocity  (cb/sec)   
  ! rvisbc......Visible beam clear sky (Refer to downward surface 
  !             shortwave fluxe)
  ! rvisdc......Visible diffuse clear sky (Refer to downward surface 
  !             shortwave fluxe)
  ! rnirbc......Near-IR beam clear skies (Refer to downward surface 
  !             shortwave fluxe)
  ! rnirdc......Near-IR diffuse clear skies (Refer to downward surface 
  !             shortwave fluxe)
  ! htrc........Heating rate due to shortwave (clear) (K/s)  
  ! dlwclr......Downward longwave (clear) At the bottom (W/m2)
  ! uswtpc......shortwave upward at top (clear) 
  ! rsclr.......net longwave at bottom (clear)                      
  ! ultclr......longwave upward at top (clear)        
  ! pie.........Constant pi=3.1415926e0   
  ! stefan......Stefan Stefan Boltzman constant    
  ! cpair.......Specific heat of air           (j/kg/k)  
  ! hl..........heat of evaporation of water     (j/kg)    
  ! grav........grav   gravity constant        (m/s**2)    
  ! snomel......Calor latente de fusao is expressed in (j m-1)
  ! tf..........Temperatura de congelamento (K)=273.16e0   
  ! clai........heat capacity of foliage	    
  ! cw..........liquid water heat capacity     (j/m**3)  
  ! gasr........gas constant of dry air        (j/kg/k)  
  ! epsfac......Constante 0.622 Razao entre as massas 
  !             moleculares do vapor e do ar seco   
  ! athird......Constant athird =                 1.0e0/3.0e0
  ! tice........tice=271.16 ice temperature ice
  ! oceald......oceald = 0.0419e0
  ! z0ice ......Roughness length of ice   
  ! icealn......near-ir beam surface albedo or near-ir diffuse surface albedo
  ! icealv......visible beam surface albedo or visible diffuse surface albedo  
  ! nfprt.......standard print out unit
  !             0 no print, 1 less detail, 2 more detail, 3 most detail
  ! nferr.......error print out unit
  !             0 no print, 1 less detail, 2 more detail, 3 most detail
  ! ifprt.......print control:
  !             0 no print, 1 less detail, 2 more detail, 3 most detail
  ! dtc3x.......time increment dt    
  ! epsflt......time filtering factor ,where {a} is time filtered value
  !             {a(t)}=a(t)+epsflt*(a(t+dt)+{a(t-dt)}-2.0*a(t))   
  ! intg........intg=1  implicit one dt integration
  !             intg=2  implicit two dt leapfrog integration    
  ! mon.........Number of month at year (1-12)     
  ! nmax........Number of point grid on continent     
  ! nsx.........Phenology dates to fall within one year period   
  ! itype.......Classe de textura do solo 
  ! vcover(iv)..Fracao de cobertura de vegetacao iv=1 Top
  ! vcover(iv)..Fracao de cobertura de vegetacao iv=2 Bottom
  ! z0x.........Roughness length   
  ! d...........Displacement height	  
  ! rdc.........constant related to aerodynamic resistance   
  ! rbc.........Constant related to bulk boundary layer resistance   
  ! z0..........Roughness length    
  ! qm..........Reference specific humidity (fourier)  
  ! tm..........Reference temperature    (fourier)                (k)  
  ! um..........Razao entre zonal pseudo-wind (fourier) e seno da
  !             colatitude  
  ! vm..........Razao entre meridional pseudo-wind (fourier) e seno da
  !             colatitude  
  ! psur........Surface pressure in mb
  ! ppc.........Precipitation rate ( cumulus )           (mm/s) 
  ! ppl.........Precipitation rate ( large scale )       (mm/s) 
  ! radn........Downward sw/lw radiation at the surface
  ! tc..........Temperatura da copa "dossel"(K)   
  ! tg .........Temperatura da superficie do solo (K)  
  ! td..........Temperatura do solo profundo (K)   
  ! capac(iv)...Agua interceptada iv=1 no dossel "water store capacity
  !             of leaves"(m)
  ! capac(iv)...Agua interceptada iv=2 na cobertura do solo (m)
  ! w(id).......Grau de saturacao de umidade do solo id=1 na camada superficial
  ! w(id).......Grau de saturacao de umidade do solo id=2 na camada de raizes
  ! w(id).......Grau de saturacao de umidade do solo id=3 na camada de drenagem
  ! ra..........Resistencia Aerodinamica (s/m)    
  ! rb..........bulk boundary layer resistance    
  ! rd..........Aerodynamic resistance between ground      (s/m)
  !             and canopy air space   
  ! rc..........Resistencia do topo da copa    
  ! rg..........Resistencia da base da copa    
  ! tcta........Diferenca entre tc-ta                      (k)  
  ! tgta........Diferenca entre tg-ta                      (k)  
  ! ta..........Temperatura no nivel de fonte de calor do dossel (K)    
  ! ea..........Pressure of vapor    
  ! etc.........Pressure of vapor at top of the copa   
  ! etg.........Pressao de vapor no base da copa   
  ! btc.........btc(i)=EXP(30.25353  -5418.0  /tc(i))/(tc(i)*tc(i)).   
  ! btg.........btg(i)=EXP(30.25353  -5418.0  /tg(i))/(tg(i)*tg(i))   
  ! u2..........wind speed at top of canopy    
  ! radt........net heat received by canopy/ground vegetation  
  ! par.........par incident on canopy   
  ! pd..........ratio of par beam to total par    
  ! rst ........Resisttencia Estomatica "Stomatal resistence" (s/m)   
  ! rsoil.......Resistencia do solo (s/m) 
  ! phroot......Soil moisture potentials in root zone of each
  !             vegetation layer and summed soil+root resistance.
  ! hrr.........rel. humidity in top layer   
  ! phsoil......soil moisture potential of the i-th soil layer
  ! cc..........heat capacity of the canopy    
  ! cg..........heat capacity of the ground    
  ! satcap......saturation liquid water capacity         (m)
  ! snow........snow amount  
  ! dtc ........dtc(i)=pblsib(i,2,5)*dtc3x  
  ! dtg.........dtg(i)=pblsib(i,1,5)*dtc3x   
  ! dtm.........dtm(i)=pblsib(i,3,5)*dtc3x   
  ! dqm ........dqm(i)=pblsib(i,4,5)*dtc3x   
  ! stm.........Variavel utilizada mo cal. da Resisttencia   
  ! extk........extinction coefficient  
  ! radfac......Fractions of downward solar radiation at surface
  !             passed from subr.radalb
  ! closs.......Radiation loss from canopy 
  ! gloss.......Radiation loss from ground 
  ! thermk......Canopy emissivity
  ! p1f   
  ! p2f   
  ! ect.........Transpiracao no topo da copa (J/m*m)   
  ! eci.........Evaporacao da agua interceptada no topo da copa (J/m*m)   
  ! egt.........Transpiracao na base da copa (J/m*m)   
  ! egi.........Evaporacao da neve (J/m*m)   
  ! egs.........Evaporacao do solo arido (J/m*m)   
  ! ec..........Soma da Transpiracao e Evaporacao da agua interceptada pelo
  !             topo da copa   ec   (i)=eci(i)+ect(i)    
  ! eg..........Soma da transpiracao na base da copa +  Evaporacao do solo arido
  !             +  Evaporacao da neve  " eg   (i)=egt(i)+egs(i)+egi(i)"    
  ! hc..........Total sensible heat lost of top from the veggies.    
  ! hg..........Total sensible heat lost of base from the veggies.    
  ! ecidif......check if interception loss term has exceeded canopy storage
  !             ecidif(i)=MAX(0.0   , eci(i)-capac(i,1)*hlat3 )
  ! egidif......check if interception loss term has exceeded canopy storage
  !             ecidif(i)=MAX(0.0   , egi(i)-capac(i,1)*hlat3 )
  ! ecmass......Mass of water lost of top from the veggies.
  ! egmass......Mass of water lost of base from the veggies.
  ! etmass......Total mass of water lost from the veggies.
  ! hflux.......Total sensible heat lost from the veggies 
  ! chf.........Heat fluxes into the canopy  in w/m**2   
  ! shf.........Heat fluxes into the ground, in w/m**2   
  ! fluxef......Modified to use force-restore heat fluxes
  !             fluxef(i) = shf(i) - cg(i)*dtg(i)*dtc3xi " Garrat pg. 227"
  ! roff........runoff (escoamente superficial e drenagem)(m)  
  ! zlwup 
  ! drag........tensao superficial  
  ! hgdtg.......n.b. fluxes expressed in joules m-2 
  ! hgdtc.......n.b. fluxes expressed in joules m-2 
  ! hgdtm.......n.b. fluxes expressed in joules m-2 
  ! hcdtg.......n.b. fluxes expressed in joules m-2 
  ! hcdtc.......n.b. fluxes expressed in joules m-2 
  ! hcdtm.......n.b. fluxes expressed in joules m-2 
  ! egdtg.......partial derivative calculation for latent heat 
  ! egdtc.......partial derivative calculation for latent heat 
  ! egdqm.......partial derivative calculation for latent heat 
  ! ecdtg.......partial derivative calculation for latent heat 
  ! ecdtc.......partial derivative calculation for latent heat 
  ! ecdqm.......partial derivative calculation for latent heat 
  ! deadtg
  ! deadtc
  ! deadqm
  ! bps   
  ! psb   
  ! dzm.........Altura media de referencia  para o vento para o calculo
  !             da estabilidade do escoamento   
  ! em..........Pressao de vapor da agua    
  ! gmt.........temperature related matrix virtual temperature tendency
  !             due to vertical diffusion   
  ! gmq.........specific humidity related matrix specific humidity of
  !             reference (fourier)   
  ! gmu.........wind related matrix   
  ! cu..........Friction  transfer coefficients.    
  ! cuni........Neutral friction transfer  coefficients.  
  ! ctni........Neutral heat transfer coefficients.  
  ! ustar.......Surface friction velocity  (m/s) 
  ! salb........surface albedo  
  ! tgeff.......effective ground temperature 
  ! cosz........Cosine of zenith angle  
  ! sinclt......sinclt=SIN(colrad(latco))"seno da colatitude"
  ! rhoair......Desnsidade do ar
  ! psy.........(cp/(hl*epsfac))*psur(i)   
  ! rcp.........densidade do ar vezes o calor especifico do ar   
  ! wc..........Minimo entre 1 e a razao entre a agua interceptada pelo
  !             indice de area foliar no topo da copa    
  ! wg..........Minimo entre 1 e a razao entre a agua interceptada pelo
  !             indice de area foliar na base da copa    
  ! fc..........Condicao de oravalho 0 ou 1 na topo da copa    
  ! fg..........Condicao de oravalho 0 ou 1 na base da copa    
  ! hrr.........rel. humidity in top layer    
  ! icheck......this version assumes dew-free conditions "icheck=1" to
  !              estimate ea for buoyancy term in vntmf or ra.
  ! ssib
  ! yvisb.......Downward Surface shortwave fluxe visible beam (cloudy)  
  ! yvisd.......Downward Surface shortwave fluxe visible diffuse (cloudy) 
  ! ynirb.......Downward Surface shortwave fluxe Near-IR beam (cloudy) 
  ! ynird.......Downward Surface shortwave fluxe Near-IR diffuse (cloudy) 
  ! yswtop......swinc....solar input at top of atmosphere
  ! yvisbc......Downward Surface shortwave fluxe visible beam (clear) 
  ! yvisdc......Downward Surface shortwave fluxe visible diffuse (clear)
  ! ynirbc......Downward Surface shortwave fluxe Near-IR beam (clear)
  ! ynirdc......Downward Surface shortwave fluxe Near-IR diffuse (clear)
  ! cldsav......Cloud cover
  ! cp..........Specific heat of air           (j/kg/k)      
  ! hl..........heat of evaporation of water     (j/kg)      
  ! rgas........gas constant of dry air        (j/kg/k)   
  ! g...........grav   gravity constant        (m/s**2) 	   
  ! solcon......solar constant (wgne value)    (w/m**2)  
  ! rmwmd.......fracao molar entre a agua e o ar   
  ! swint.......sw subr. call interval in hours
  !             swint has to be less than or equal to trint
  !                              and mod(trint,swint)=0
      
  ! trint.......ir subr. call interval in hours      
  ! yrl.........length of year in days        
  ! idate(4)....output : idate(1) = initial hour of
  !                      idate(2) = day of month.
  !                      idate(3) = month of year.
  !                      idate(4) = year.
 
  ! idatec(4)...output : idatec(1)=current hour of
  !			 idatec(2)=current day of month.
  !			 idatec(3)=current month of year.
  !			 idatec(4)=current year.
 
  ! kt..........hour of present  time step       
  ! ktm.........hour of previous time step       
  ! jdt.........time step in getdia      
  ! monl(12)....length of each month in days   
  ! irad........the physical process radiation 	 
  !                     irad = yes   included
  !                     irad = no    excluded
  ! iccon.......the physical process cumulus convection(kuo)
  !                     iccon = yes   included
  !                     iccon = no    excluded	     
  ! icld........>>> icld = 1    : old cloud emisivity (optical depth) setting      
  !		       ccu :	 0.05 *dp				     
  !		       css :	 0.025*dp	      for ice cloud t<253.0  
  !				 0.05 *dp	      for ice cloud t>253.0  
  !		>>> icld = 2    : new cloud emisivity (optical depth) setting      
  !		       ccu :	 (0.16)*dp				     
  !		       css :	  0.0			      t<-82.5c       
  !				 (2.0e-6*(t-tcrit)**2)*dp    -82.5<t<-10.0c  
  !				 (6.949e-3*(t-273)+.08)*dp   -10.0<t< 0.0c   
  !				 (0.08)*dp		     -10.0<t< 0.0c   
  !		>>> icld = 3    : ccm3 based cloud emisivity 		     
  ! inalb.......inalb Input two types of surface albedo 
  !                >>> inalb= 1 : input two  types surfc albedo (2 diffused)
  !                             direct beam albedos are calculated by the subr.
  !                >>> inalb= 2 : input four types surfc albedo (2 diff,2 direct)
	     
  ! co2val......co2val is wgne standard value in ppm "co2val = /345.0/	     
  ! delt........time interval in sec (fixed throuh the integration)	     
  ! filta.......weight used on central time
  !              step of robert time filter.
  !              set in main routine "smf".filta=0.92e0	     
  ! nfin0.......input  file at time level t-dt	     
  ! nfin1.......input  file at time level t     
  ! initlz......constant initlz=2.	     
  ! nfcnv0......initial information on convective clouds for int. radiation	     
  ! nfcldr......constant nfcldr = 74    
  ! tbase.......from common commnt constant tbase =  273.15e00	  
  ! latco.......latitude 	  
  ! dodia.......Variable logical for search for combined field components. 
  ! lvavl.......levels in available diagnostic (1 or kmax)
  ! nuavl.......unit code of available diagnostic
  ! itavl.......type of available diagnostic (1 gaussian, 2 spectral)
  ! iavrq.......Number of requested diagnostic	
  ! ixavl.......Number available diagnostic components for combined fields
  ! inavl.......Number available diagnostic similar requested diagnostic
  ! iclcd.......requested diagnostic calculation code (0 direct
  !             calculation, > 0 add to requested field number iclcd,
  !             < 0 subtract from requested field number -iclcd )
  ! nucf........nurq  = unit code of requested diagnostic	 
  ! ixcf........Number of requested diagnostic	   
  ! incf........combined fields 
  ! kravl.......Number of available diagnostics equivalence 
  !              the desired diagnostic
  ! krcf........combined fields 
  ! jrcf........combined fields 
  ! ndt.........time step in getdia 	     
  ! latitu......latco grid point reference the latitude 	     
  ! cdhl........logical indicator for dhn output prognostics
  ! ustr........surface zonal stress umom(i)=fmom*um(ncount),
  !               where .fmom  momentum flux      in n/m**2
  !               fmom= rhoair(ncount)*cu(ncount)*ustar(ncount)
  !               um  (ncount)=gu (i,1)/sinclt
  !               gu          = (zonal velocity)*sin(colat)

  ! vstr........surface meridional stress.vmom(i)=rho(i)*gv(i)*rmi(i)
  !                           rho  (i)=gps(i)/(gr100*gt(i))
  !                           gr100 =gasr*0.01
  ! first.......control logical variable .true. or .false. 
  ! mxrdcc......use maximum random converage for radiative conv. clouds
  !               constant logical mxrdcc = .true.
  ! lcnvl.......the lowest layer index where non-convective clouds can
  !               occur (ben says this should be 2 or more)
  !               constant lcnvl = 2 
  ! lthncl......Minimum depth in mb of non-zero low level cloud
  !             consta lthncl=80
  ! convc.......ncols*jmax convective cloud cover in 3 hr. avrage
  ! convt.......ncols*jmax convective cloud top  (sigma layer) 
  ! convb.......ncols*jmax convective cloud base (sigma layer)
  ! convts
  ! convcs
  ! convbs
  ! sigki ......sigki (k)=1.0e0/EXP(akappa*LOG(sig(k))),  where "sig"
  !             sigma coordinate at middle of layer and akappa=gasr/cp 
  ! xvisb.......Downward Surface shortwave fluxe visible beam (cloudy)  
  ! xvisd.......Downward Surface shortwave fluxe visible diffuse (cloudy) 
  ! xnirb.......Downward Surface shortwave fluxe Near-IR beam (cloudy)  
  ! xnird.......Downward Surface shortwave fluxe Near-IR diffuse (cloudy)   
  ! xswtop......shortwave upward at top  or  shortwave upward at top (clear)
  ! xvisbc......Downward Surface shortwave fluxe visible beam (clear) 
  ! xvisdc......Downward Surface shortwave fluxe visible diffuse (clear)
  ! xnirbc......Downward Surface shortwave fluxe Near-IR beam (clear)
  ! xnirdc......Downward Surface shortwave fluxe Near-IR diffuse (clear)
  !==========================================================================
  !
  
  IMPLICIT NONE
  INTEGER, INTENT(in   ) :: ncols
  INTEGER, INTENT(in   ) :: kmax
  REAL   , PARAMETER :: alon   =  0.0 
  INTEGER, PARAMETER :: ityp   = 13
  INTEGER, PARAMETER :: imon   = 12
  INTEGER, PARAMETER :: icg    =  2
  INTEGER, PARAMETER :: iwv    =  3
  INTEGER, PARAMETER :: ild    =  2
  INTEGER, PARAMETER :: idp    =  3
  INTEGER, PARAMETER :: ibd    =  2

    
  REAL,    INTENT(inout) :: gt    (ncols,kmax)
  REAL,    INTENT(inout) :: gq    (ncols,kmax)
  REAL,    INTENT(in   ) :: gu    (ncols,kmax)
  REAL,    INTENT(in   ) :: gv    (ncols,kmax)
  REAL,    INTENT(in   ) :: gps   (ncols)      
  REAL,    INTENT(inout) :: tc0   (ncols)
  REAL,    INTENT(inout) :: tg0   (ncols)
  REAL,    INTENT(inout) :: td0   (ncols)
  REAL,    INTENT(inout) :: w0    (ncols,3)
  REAL,    INTENT(inout) :: capac0(ncols,2)
  REAL,    INTENT(inout) :: tcm   (ncols)
  REAL,    INTENT(inout) :: tgm   (ncols)
  REAL,    INTENT(inout) :: tdm   (ncols)
  REAL,    INTENT(inout) :: wm    (ncols,3)
  REAL,    INTENT(inout) :: capacm(ncols,2)
  REAL,    INTENT(in   ) :: ppli  (ncols)
  REAL,    INTENT(in   ) :: ppci  (ncols)
  REAL,    INTENT(inout) :: gyu   (ncols,kmax)
  REAL,    INTENT(inout) :: gyv   (ncols,kmax)
  REAL,    INTENT(inout) :: gtd   (ncols,kmax)
  REAL,    INTENT(inout) :: gqd   (ncols,kmax)
  REAL,    INTENT(out  ) :: htr   (ncols,kmax)
  REAL,    INTENT(out  ) :: clr   (ncols,kmax)
  INTEGER, INTENT(in   ) :: imask (ncols)
  REAL,    INTENT(out  ) :: rvisb (ncols)
  REAL,    INTENT(out  ) :: rvisd (ncols)
  REAL,    INTENT(out  ) :: rnirb (ncols)
  REAL,    INTENT(out  ) :: rnird (ncols)
  REAL,    INTENT(out  ) :: dlwbot(ncols)
  REAL,    INTENT(inout) :: tsea  (ncols)
  REAL,    INTENT(in   ) :: colrad(ncols)
  REAL,    INTENT(in   ) :: lonrad(ncols)
  REAL,    INTENT(in   ) :: sig   (kmax)
  REAL,    INTENT(in   ) :: sigml (kmax)
  REAL,    INTENT(in   ) :: delsig(kmax)
  INTEGER, INTENT(in   ) :: ifday
  REAL,    INTENT(in   ) :: tod   
  REAL,    INTENT(out  ) :: avisb (ncols)
  REAL,    INTENT(inout) :: avisd (ncols)
  REAL,    INTENT(out  ) :: anirb (ncols)
  REAL,    INTENT(out  ) :: anird (ncols)
  REAL,    INTENT(out  ) :: dswtop(ncols)
  REAL,    INTENT(out  ) :: rs    (ncols)
  REAL,    INTENT(out  ) :: ulwtop(ncols)
  REAL,    INTENT(inout) :: gl0   (ncols)
  REAL,    INTENT(inout) :: zorl  (ncols)
  REAL,    INTENT(in   ) :: gtt   (ncols,kmax)
  REAL,    INTENT(in   ) :: gqq   (ncols,kmax)
  REAL,    INTENT(in   ) :: sheleg(ncols)
  REAL,    INTENT(inout) :: tseam (ncols)
  REAL,    INTENT(in   ) :: omg   (ncols,kmax)
  REAL,    INTENT(out  ) :: rvisbc(ncols)
  REAL,    INTENT(out  ) :: rvisdc(ncols)
  REAL,    INTENT(out  ) :: rnirbc(ncols)
  REAL,    INTENT(out  ) :: rnirdc(ncols)
  REAL,    INTENT(out  ) :: htrc  (ncols,kmax)
  REAL,    INTENT(out  ) :: dlwclr(ncols)
  REAL,    INTENT(out  ) :: uswtpc(ncols)
  REAL,    INTENT(out  ) :: rsclr (ncols)
  REAL,    INTENT(out  ) :: ultclr(ncols)

  INTEGER :: mon		      ! from common sibtim
  INTEGER :: nmax		      ! from common sibtim

  REAL,    INTENT(in   ) :: ssib(ncols)            ! from common comsib
  !
  !     this is for interpolating shortwave rad at ground
  !
  REAL,    INTENT(out  ) :: yvisb (ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: yvisd (ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: ynirb (ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: ynird (ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: yswtop(ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: yvisbc(ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: yvisdc(ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: ynirbc(ncols)	 ! from common comswint
  REAL,    INTENT(out  ) :: ynirdc(ncols)	 ! from common comswint
  REAL,    INTENT(inout) :: cldsav(ncols)	 ! from common comswint
  !
  !     these are for monitoring of gpv in gfidi.
  !

  INTEGER, INTENT(in   ) :: latco                      ! from common commnt

  REAL,    INTENT(out  ) :: ustr(ncols)            ! from common cmndhn 
  REAL,    INTENT(out  ) :: vstr(ncols)            ! from common cmndhn 

  REAL,    INTENT(in   ) :: convc (ncols)	       ! from common  cumcnv
  REAL,    INTENT(in   ) :: convt (ncols)	       ! from common  cumcnv
  REAL,    INTENT(in   ) :: convb (ncols)	       ! from common  cumcnv

  REAL,    INTENT(in   ) :: convts(ncols)          ! from common shalcom
  REAL,    INTENT(in   ) :: convcs(ncols)          ! from common shalcom
  REAL,    INTENT(in   ) :: convbs(ncols)          ! from common shalcom

  REAL,    INTENT(out) :: sigki (kmax)
  REAL,    INTENT(out) :: xvisb (ncols)
  REAL,    INTENT(out) :: xvisd (ncols)
  REAL,    INTENT(out) :: xnirb (ncols)
  REAL,    INTENT(out) :: xnird (ncols)
  REAL,    INTENT(out) :: xswtop(ncols)
  REAL,    INTENT(out) :: xvisbc(ncols)
  REAL,    INTENT(out) :: xvisdc(ncols)
  REAL,    INTENT(out) :: xnirbc(ncols)
  REAL,    INTENT(out) :: xnirdc(ncols)  

  REAL,    INTENT(inout) :: ps    (ncols)
  REAL   , INTENT(inout) :: var   (ncols)               ! from common sf
  REAL,    INTENT(out  ) :: sens  (ncols)
  REAL,    INTENT(out  ) :: evap  (ncols)    
  REAL   , INTENT(IN   ) :: cos2d  (ncols)    
  LOGICAL, INTENT(IN   ) :: intcosz
  INTEGER :: imx
  INTEGER :: jmax
  INTEGER :: icld			! from common comctm
  INTEGER :: inalb			! from common comctm
  REAL                   :: uswtop(ncols)



  INTEGER :: nsx			! from common nwhea9

  INTEGER :: itype(ncols)		 ! from common comtyp
  REAL :: zlt2    (ncols,icg)
  REAL :: green2  (ncols,icg)
  REAL :: chil2   (ncols,icg)
  REAL :: vcover  (ncols,icg)
  REAL :: topt2   (ncols,icg)
  REAL :: tll2	  (ncols,icg)
  REAL :: tu2	  (ncols,icg)
  REAL :: defac2  (ncols,icg)
  REAL :: ph12	  (ncols,icg)
  REAL :: ph22	  (ncols,icg)
  REAL :: rstpar2(ncols,icg,iwv)        ! from common comres
  REAL :: z0x	(ncols) 	      ! from common comres
  REAL :: d	(ncols) 	      ! from common comres
  REAL :: rdc	(ncols) 	      ! from common comres
  REAL :: rbc	(ncols) 	      ! from common comres
  REAL :: z0	(ncols) 	      ! from common comres
  !
  !     the size of working area is ncols*187
  !     atmospheric parameters as boudary values for sib
  !
  REAL :: qm  (ncols)		      ! from common comatm
  REAL :: tm  (ncols)		      ! from common comatm
  REAL :: um  (ncols)		      ! from common comatm
  REAL :: vm  (ncols)		      ! from common comatm
  REAL :: psur(ncols)		      ! from common comatm
  REAL :: ppc (ncols)		      ! from common comatm
  REAL :: ppl (ncols)		      ! from common comatm
  REAL :: radn(ncols,3,2)	      ! from common comatm
  !
  !     prognostic variables
  !
  REAL :: tc   (ncols)  	      ! from common comprg 
  REAL :: tg   (ncols)  	      ! from common comprg 
  REAL :: td   (ncols)  	      ! from common comprg 
  REAL :: capac(ncols,2)	      ! from common comprg 
  REAL :: w    (ncols,3)	      ! from common comprg 
  !
  !     variables calculated from above and ambient conditions
  !
  REAL :: ra	(ncols) 	    ! from common comdia
  REAL :: rb	(ncols) 	    ! from common comdia
  REAL :: rd	(ncols) 	    ! from common comdia
  REAL :: rc	(ncols) 	    ! from common comdia
  REAL :: rg	(ncols) 	    ! from common comdia
  REAL :: tcta  (ncols) 	    ! from common comdia
  REAL :: tgta  (ncols) 	    ! from common comdia
  REAL :: ta	(ncols) 	    ! from common comdia
  REAL :: ea	(ncols) 	    ! from common comdia
  REAL :: etc	(ncols) 	    ! from common comdia
  REAL :: etg	(ncols) 	    ! from common comdia
  REAL :: btc	(ncols) 	    ! from common comdia
  REAL :: btg	(ncols) 	    ! from common comdia
  REAL :: u2	(ncols) 	    ! from common comdia
  REAL :: radt  (ncols,icg)	    ! from common comdia
  REAL :: par	(ncols,icg)	    ! from common comdia
  REAL :: pd	(ncols,icg)	    ! from common comdia
  REAL :: rst	(ncols,icg)	    ! from common comdia
  REAL :: rsoil (ncols) 	    ! from common comdia
  REAL :: phroot(ncols,icg)	    ! from common comdia
  REAL :: hrr	(ncols) 	    ! from common comdia
  REAL :: phsoil(ncols,idp)	    ! from common comdia
  REAL :: cc	(ncols) 	    ! from common comdia
  REAL :: cg	(ncols) 	    ! from common comdia
  REAL :: satcap(ncols,icg)	    ! from common comdia
  REAL :: snow  (ncols,icg)	    ! from common comdia
  REAL :: dtc	(ncols) 	    ! from common comdia
  REAL :: dtg	(ncols) 	    ! from common comdia
  REAL :: dtm	(ncols) 	    ! from common comdia
  REAL :: dqm	(ncols) 	    ! from common comdia
  REAL :: stm	(ncols,icg)	    ! from common comdia
  REAL :: extk  (ncols,icg,iwv,ibd)  ! from common comdia
  REAL :: radfac(ncols,icg,iwv,ibd)  ! from common comdia
  REAL :: closs (ncols) 	    ! from common comdia
  REAL :: gloss (ncols) 	    ! from common comdia
  REAL :: thermk(ncols) 	    ! from common comdia
  REAL :: p1f	(ncols) 	    ! from common comdia
  REAL :: p2f	(ncols) 	    ! from common comdia
  !
  !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
  !
  REAL :: ect	(ncols) 	      ! from common comflx
  REAL :: eci	(ncols) 	      ! from common comflx
  REAL :: egt	(ncols) 	      ! from common comflx
  REAL :: egi	(ncols) 	      ! from common comflx
  REAL :: egs	(ncols) 	      ! from common comflx
  REAL :: ec	(ncols) 	      ! from common comflx
  REAL :: eg	(ncols) 	      ! from common comflx
  REAL :: hc	(ncols) 	      ! from common comflx
  REAL :: hg	(ncols) 	      ! from common comflx
  REAL :: ecidif(ncols) 	      ! from common comflx
  REAL :: egidif(ncols) 	      ! from common comflx
  REAL :: ecmass(ncols) 	      ! from common comflx
  REAL :: egmass(ncols) 	      ! from common comflx
  REAL :: etmass(ncols) 	      ! from common comflx
  REAL :: hflux (ncols) 	      ! from common comflx
  REAL :: chf	(ncols) 	      ! from common comflx
  REAL :: shf	(ncols) 	      ! from common comflx
  REAL :: fluxef(ncols) 	      ! from common comflx
  REAL :: roff  (ncols) 	      ! from common comflx
  REAL :: zlwup (ncols) 	      ! from common comflx
  REAL :: drag  (ncols) 	      ! from common comflx
  !
  !     this is for coupling with closure turbulence model
  !
  REAL :: bps	(ncols) 	      ! from common comcup  
  REAL :: psb	(ncols) 	      ! from common comcup  
  REAL :: dzm	(ncols) 	      ! from common comcup  
  REAL :: em	(ncols) 	      ! from common comcup  
  REAL :: gmt	(ncols,3)	      ! from common comcup  
  REAL :: gmq	(ncols,3)	      ! from common comcup  
  REAL :: gmu	(ncols,4)	      ! from common comcup  
  REAL :: cu	(ncols) 	      ! from common comcup  
  REAL :: cuni  (ncols) 	      ! from common comcup  
  REAL :: ctni  (ncols) 	      ! from common comcup  
  REAL :: ustar (ncols) 	      ! from common comcup  
  REAL :: salb  (ncols,2,2)	      ! from common comcup  
  REAL :: tgeff (ncols) 	      ! from common comcup  
  REAL :: cosz  (ncols) 	      ! from common comcup  
  REAL :: sinclt		      ! from common comcup  
  REAL :: rhoair(ncols) 	      ! from common comcup 
  REAL :: psy	(ncols) 	      ! from common comcup  
  REAL :: rcp	(ncols) 	      ! from common comcup  
  REAL :: wc	(ncols) 	      ! from common comcup  
  REAL :: wg	(ncols) 	      ! from common comcup  
  REAL :: fc	(ncols) 	      ! from common comcup  
  REAL :: fg	(ncols) 	      ! from common comcup  
  REAL :: hr	(ncols) 	      ! from common comcup 

  INTEGER :: icheck(ncols)            ! from common comchk
  
  INTEGER :: ndt                      ! from common tempor 
  INTEGER :: latitu                   ! from common tempor  
  
  REAL    :: swrk  (ncols,57)                            ! from common sibwrk
  REAL    :: tmtx  (ncols,kmax,3)
  REAL    :: qmtx  (ncols,kmax,3)
  REAL    :: umtx  (ncols,kmax,4)
  REAL    :: bfrg  (ncols)
  REAL    :: bfr1  (ncols)
  REAL    :: bfr2  (ncols)
  REAL    :: gtmp  (ncols,kmax)   
  REAL    :: gwv   (ncols,kmax)
  REAL    :: rnet  (ncols)
  REAL    :: q     (ncols,kmax)   
  REAL    :: pcb   (ncols)      
  REAL    :: xsea  (ncols)
  REAL    :: relhum(ncols,kmax)
  REAL    :: tsurf (ncols)
  REAL    :: qsurf (ncols)
  REAL    :: rmi   (ncols)
  REAL    :: rhi   (ncols)
  REAL    :: umom  (ncols)
  REAL    :: vmom  (ncols)
  REAL    :: strcut
  REAL    :: slrad (ncols)
  REAL    :: zenith(ncols)
  REAL    :: cond  (ncols)
  REAL    :: stor  (ncols)
  
  REAL    :: tmin  (ncols)
  REAL    :: tmax  (ncols)
  REAL    :: cldtot(ncols,kmax)
  REAL    :: cldinv(ncols,kmax)
  REAL    :: cldsat(ncols,kmax)
  REAL    :: cldcon(ncols,kmax)
  REAL    :: cldson(ncols,kmax)
  REAL    :: clwd  (ncols,kmax)
  REAL    :: emisd (ncols,kmax)
  REAL    :: taud  (ncols,kmax)
  REAL    :: work  (ncols)
  !
  !     the names of some variables in comphc have to be changed to
  !     avoid conflict with names in sibcon (comsibw)
  !
  !     intg=1  implicit one dt integration
  !     =2  implicit two dt leapfrog integration
  !     epsflt  time filtering factor ,where {a} is time filtered value
  !     {a(t)}=a(t)+epsflt*(a(t+dt)+{a(t-dt)}-2.0*a(t))
  !     
  LOGICAL, SAVE :: initl=.TRUE.
  LOGICAL, SAVE :: cloudr
  REAL,    SAVE :: pai12
  REAL,    SAVE :: akappa
  REAL,    SAVE :: sdelt
  REAL,    SAVE :: ratio
  REAL,    SAVE :: etime
  REAL,    SAVE :: xday

  LOGICAL :: ghl_local
  INTEGER :: k 
  INTEGER :: i 
  INTEGER :: j 
  INTEGER*4 :: lrecl 
  INTEGER :: ncount 
  INTEGER :: ntyp 
  INTEGER :: jlat 
  INTEGER :: ind 
  INTEGER :: itr 
  INTEGER :: intstp 
  INTEGER :: modstp 
  INTEGER :: mmax
  REAL    :: csqiv2 (ncols)
  REAL    :: deltm 
  REAL    :: sindel 
  REAL    :: cosdel 
  REAL    :: fimxi
  REAL    :: ctime 
  REAL    :: cos2  (ncols)
  REAL    :: frh 
  REAL    :: btime 
  REAL    :: atime 
  REAL    :: s0 
  REAL    :: cpsy 
  REAL    :: rbyg 
  REAL    :: r100 
  REAL    :: f
  REAL    :: ocealb 
  REAL    :: tice01 
  REAL    :: radint 
  REAL    :: fstp 
  REAL    :: fstp1 
  REAL    :: fmom 
  REAL    :: rfac 
  REAL    :: dtc3xi
  REAL    :: one   =1.0e0 
  REAL    :: thousd=1.0e3
  REAL    :: f3600 =3.6e3                    
  REAL    :: cnvfac=14.3353e-04
  REAL    :: tenth =0.1e0
  REAL    :: delq  =0.608e0
  INTEGER :: month   (ncols) 
  INTEGER :: month2  (ncols) 
  INTEGER :: nn 
  REAL    :: colrad2 (ncols)
  REAL    :: zenith1 (ncols) 
  REAL    :: zenith2 (ncols)
  REAL    :: sinclt2 (ncols)
  REAL    :: cos3 
  REAL    :: chug  (ncols,kmax)
  REAL    :: chvg  (ncols,kmax)
  REAL    :: xdrag (ncols)
  REAL    :: ydrag (ncols)
  REAL    :: julday
!--------------------------------------------------------------------------
!	 locations for available diagnostics in this subroutine
!--------------------------------------------------------------------------     
!  INTEGER :: ndshf  = 16   ! sensible heat flux 		    
!  INTEGER :: ndlhf  = 17   ! latent heat flux				    
!  INTEGER :: ndustr = 18   ! surface zonal stress		    
!  INTEGER :: ndvstr = 19   ! surface meridional stress  	    
!  INTEGER :: nddtvd = 39   ! vertical diffusion heating 	    
!  INTEGER :: ndduvd = 41   ! vertical diffusion zonal momentum change	    
!  INTEGER :: nddvvd = 42   ! vertical diffusion meridional momentum change
!  INTEGER :: nddqvd = 40   ! vertical diffusion moistening	    
!  INTEGER :: ndtsb  =  6   ! time mean surface temperature	    
!  INTEGER :: ndcld  = 20   ! cloud cover			    
!  INTEGER :: nddtlw = 31   ! longwave heating				    
!  INTEGER :: nddtsw = 32   ! shortwave heating  		    
!  INTEGER :: ndlwdb = 21   ! longwave downward at bottom	    
!  INTEGER :: ndlwub = 22   ! longwave upward at bottom  	    
!  INTEGER :: ndlwut = 23   ! longwave upward at top		    
!  INTEGER :: ndswdt = 24   ! shortwave downward at top  	    
!  INTEGER :: ndswdb = 25   ! shortwave downward at ground	    
!  INTEGER :: ndswea = 28   ! shortwave absorbed by the earth/atmosphere 
!  INTEGER :: ndswe  = 29   ! shortwave absorbed by the ground		    
!  INTEGER :: ndswut = 27   ! shortwave upward at top			    
!  INTEGER :: ndswub = 26   ! shortwave upward at bottom 	    
!  INTEGER :: ndlwnb = 30   ! net longwave at bottom		    
!  INTEGER :: ndintl = 15   ! interception loss  		    
!  INTEGER :: ndrnof = 13   ! runoff				    
!  INTEGER :: ndldbc = 56   ! longwave downward at bottom (clear)    
!  INTEGER :: ndlutc = 57   ! longwave upward at top (clear)	    
!  INTEGER :: ndsdbc = 58   ! shortwave downward at ground (clear)   
!  INTEGER :: ndseac = 61   ! shortwave absorbed by the earth/atmosphere (clear)
!  INTEGER :: ndswec = 62   ! shortwave absorbed by the ground (clear)	    
!  INTEGER :: ndsutc = 60   ! shortwave upward at top (clear)		    
!  INTEGER :: ndsubc = 59   ! shortwave upward at bottom (clear)     
!  INTEGER :: ndlnbc = 63   ! net longwave at bottom (clear)	    
!  INTEGER :: ndta   = 67   ! temperature of canopy air space		    
!  INTEGER :: ndea   = 68   ! vapor pressure of canopy air space     
!  INTEGER :: ndlhbs = 69   ! bare soil latent heat		    
!  INTEGER :: ndctd  = 72   ! vertical dist total cloud cover		    
!  INTEGER :: ndcin  = 73   ! inversion cloud				    
!  INTEGER :: ndcst  = 74   ! supersaturation cloud		    
!  INTEGER :: ndccn  = 75   ! convective cloud				    
!  INTEGER :: ndcsn  = 76   ! shallow convective cloud			    
!  INTEGER :: ndclw  = 77   ! cloud liquid water path			    
!  INTEGER :: ndemi  = 78   ! longwave cloud emissivity  	    
!  INTEGER :: ndsod  = 79   ! shortwave cloud optical depth	    
!  INTEGER :: ndrm   = 80   ! momentum flux resistance			    
!  INTEGER :: ndrarh = 81   ! canopy air spc to ref. lvl resistance  
!  INTEGER :: ndrb   = 82   ! canopy air spc to canopy resistance    
!  INTEGER :: ndrd   = 83   ! canopy air spc to ground resistance    
!  INTEGER :: ndrg   = 84   ! ground cover resistance			    
!  INTEGER :: ndrs   = 85   ! bare soil surface resistance	    
!  INTEGER :: ndtxgw = 43   ! gravity wave drag surface zonal stress 
!  INTEGER :: ndtygw = 44   ! gravity wave drag surface meridional stress
!  INTEGER :: ndugwd = 45   ! gravity wave drag zonal momentum change	    
!  INTEGER :: ndvgwd = 46   ! gravity wave drag meridional momentum change
!--------------------------------------------------------------------------
   
  ghl_local = IsGridHistoryOn() 
  if(dogwd.eq.0)then
      !
      !     change on gwdd by cptec on 29 july 1994 to improve
      !      vectorization performance on gwdd
      !
       call InitGwdd(sigml, kmax, nfprt, nferr)
       call Gwdd(ps ,gu   ,gv   ,gt   , chug, chvg, xdrag, ydrag, &
                 nfprt, var, varcut, sigml, sig, delsig, ncols, kmax)
        do  k=1, kmax
          do i=1, ncols	  	  
          ! gyu(i,k)=gyu(i,k)+chvg(i,k)
          ! gyv(i,k)=gyv(i,k)-chug(i,k)
           gyu(i,k) = gyu(i,k) - chug(i,k)
           gyv(i,k) = gyv(i,k) - chvg(i,k)
          end do
        end do
  endif
  vcover  = 0.0
  zlt2    = 0.0
  green2  = 0.0
  chil2   = 0.0
  rstpar2 = 0.0
  topt2   = 0.0
  tll2    = 0.0
  tu2	  = 0.0
  defac2  = 0.0
  ph12    = 0.0
  ph22    = 0.0
   
  if(dodia(ndtxgw))call updia(xdrag,ndtxgw,latco)
  if(dodia(ndtygw))call updia(ydrag,ndtygw,latco)
  if(dodia(ndugwd))call updia(chug,ndugwd,latco)
  if(dodia(ndvgwd))call updia(chvg,ndvgwd,latco)
   ndt=jdt
   latitu=latco
  IF(initl) THEN
      initl=.FALSE.
      intg =intgr
      pai12=pie/12.0e0
      epsflt=0.5e0 *(1.0e0 -filta)
      IF(intg.EQ.1) THEN
        epsflt=0.0e0
      END IF
      akappa=gasr/cp
      DO k=1,kmax
        sigki (k)=1.0e0/EXP(akappa*LOG(sig(k)))
      END DO
      DO i=1,ncols
        xvisd (i)=0.0e0
        xvisb (i)=0.0e0
        xnird (i)=0.0e0
        xnirb (i)=0.0e0
        xvisdc(i)=0.0e0
        xvisbc(i)=0.0e0
        xnirdc(i)=0.0e0
        xnirbc(i)=0.0e0
        xswtop(i)=0.0e0
      END DO 
      cloudr=dodia(ndctd).OR.dodia(ndcin).OR.dodia(ndcst).OR. &
             dodia(ndccn).OR.dodia(ndcsn).OR.dodia(ndclw).OR. &
             dodia(ndemi).OR.dodia(ndsod).OR.(ghl_local.AND.dogrh(ngctd,latco))
      IF(cloudr.AND.jdt.EQ.1.AND.ktm.EQ.-1.AND.kt.EQ.0)THEN
        INQUIRE (IOLENGTH=lrecl) cldtot,cldinv,cldsat,cldcon, &
             cldson,clwd,emisd,taud
        OPEN(nfcldr,access='DIRECT',form='UNFORMATTED',recl=lrecl, &
             status='UNKNOWN')
      END IF
  END IF
  IF(cloudr)THEN
      DO k=1,kmax
        DO i=1,ncols
           cldtot(i,k)=0.0e0
           cldinv(i,k)=0.0e0
           cldsat(i,k)=0.0e0
           cldcon(i,k)=0.0e0
           cldson(i,k)=0.0e0
           clwd(i,k)=0.0e0
           emisd(i,k)=0.0e0
           taud(i,k)=0.0e0
        END DO
      END DO
  ENDIF
  
  dtc3x   =dt*float(intg)

  DO i=1,ncols,1 
      csqiv2(i)   = 1.0e0/sin(colrad(i))**2
  END DO
 !
 !
 !
  DO k=1,kmax
      DO i=1,ncols
        gq(i,k)=MAX(1.0e-12,gq(i,k))
        gt(i,k)=gt(i,k)/(1.0e0+0.608e0*gq(i,k))
      END DO
  END DO
  DO k=1,kmax
      DO i=1,ncols
        gtmp(i,k)=gtt(i,k)
      END DO
  END DO  
  DO k=1,kmax
      DO i=1,ncols
        gwv (i,k)=gqq(i,k)
      END DO
  END DO
  DO i=1,ncols
      pcb(i)=tenth*gps(i)
  END DO
  IF(initlz.GE.0.AND.kt.EQ.0.AND.jdt.EQ.1) THEN
      ncount=0
      DO i=1,ncols
        IF(imask(i).GE.1) THEN
           ncount=ncount+1
           tc0(i)=gt(i,1)
           tg0(i)=gt(i,1)
           tcm(i)=gt(i,1)
           tgm(i)=gt(i,1)
           IF(sheleg(i).GT.0.0e0) THEN
              tg0(i)=MIN(tg0(i),tf-0.01e0)
              tgm(i)=MIN(tgm(i),tf-0.01e0)
           END IF
        END IF
      END DO
  END IF
  ncount=0
  DO i=1,ncols
      rhi(i)=1.0e-6
      rmi(i)=1.0e-6
      IF(imask(i).GE.1) THEN
        ncount=ncount+1
        itype(i)=imask(i)
      END IF
  END DO
  nmax=ncount
  !
  !     yamada-mellor pbl parameterization
  !
  deltm =0.5e0 *dtc3x
  CALL ympbl0( &
              gu    ,gv    ,gt    ,gq    ,delsig,ncols, &
	      kmax  ,deltm ,csqiv2,tmtx  ,qmtx  ,umtx  ,gl0   )
  !
  !     computation of astronomical parameters
  !     sdelt ;solar inclination
  !     etime ;correction factor to local time
  !     ratio ;factor relating to the distance between the earth and the sun
  !
  CALL radtim(idatec, &! INTENT(in ) :: id(4)
              sdelt , &! INTENT(out) :: delta
	      ratio , &! INTENT(out) :: ratio
	      etime , &! INTENT(out) :: etime
	      tod   , &! INTENT(in ) :: tod
	      xday  , &! INTENT(out) :: xday
	      nfprt , &! INTENT(in ) :: nfprt
              nferr , &! INTENT(in ) :: nferr
	      yrl     )! INTENT(in ) :: yrl
  !     
  !     mon is the month used for vegetation data input
  ! 
  DO i=1,ncols,1
   month(i)=idatec(2)
   IF((((colrad(i)*180.0)/3.1415926e0)-90.0)  > 0.0 ) THEN
       month(i)  =  month(i) + 6
       IF(month(i).GE.13) month(i) = month(i)-12
   END IF
  END DO
  
  ncount=0
  
  DO i=1,ncols
      IF(imask(i).GE.1) THEN
        ncount=ncount+1
        month2(i)=month(i)
	colrad2(i)=colrad(i)
	sinclt2(i) = sin(colrad(i))
      END IF
  END DO  
  
  sindel=SIN(sdelt)
  cosdel=COS(sdelt)
  fimxi  =24.0e0 /360.0
  ctime  =alon/15.0e0
  cos2   =0.0e0
  cos3   =0.0e0
  ncount =0
  frh=( MOD(tod+0.03125,f3600)-0.03125)/f3600

  DO i=1,ncols
   zenith1(i)  = sindel*cos(colrad(i))
  ENDDO

  DO i=1,ncols
      btime       = fimxi*lonrad(i)+ctime
      atime       = etime+pai12*(12.0-idatec(1)-frh-btime)
      zenith2(i)  = cosdel*sin(colrad(i))*COS(atime)
      zenith(i)   = zenith1(i) + zenith2(i)
      IF(zenith(i).GT.0.0e0) THEN
        ncount   =ncount+1
        !cos2     =cos2+zenith(i) !!!!mudanca forcada 
      END IF
  END DO 

  IF(ncount.EQ.0) ncount=1
  IF(intcosz)THEN
   !cos2=cos2/float(ncount)!!!!mudanca forcada 
    cos2=cos2d
  ELSE  
    cos2=zenith
  END IF    
  s0     =ratio*solcon
  !
  !     sib setting **
  !
  IF(nmax.GE.1) THEN
      DO i=1,nmax
        ntyp=itype(i)
	topt2	(i,1) =   topt  (ntyp,1)
	tll2	(i,1) =   tll   (ntyp,1)
	tu2	(i,1) =   tu    (ntyp,1)
	defac2  (i,1) =   defac (ntyp,1)
	ph12	(i,1) =   ph1   (ntyp,1)
	ph22	(i,1) =   ph2   (ntyp,1)
	DO  nn = 1, 3
          rstpar2 (i,1,nn)= rstpar(ntyp,1,nn)
        END DO
	zlt2   (i,1) =  zlt   (ntyp,month2(i),1)
	green2 (i,1) =  green (ntyp,month2(i),1)
	chil2  (i,1) =  chil  (ntyp,2)
        vcover (i,1) =  xcover(ntyp,month2(i),1)
      END DO
      xday=julday(idatec,tod,monl)
      CALL wheat (itype  (1:nmax)  ,&
                  nmax             ,&
		  colrad2(1:nmax)  ,&
		  month2 (1:nmax)  ,&
		  xday             ,&
		  yrl              ,&
                  idatec           ,&
		  monl             ,&
		  nsx              ,&
		  nmax             ,&
		  vcover (1:nmax,:),&
		  zlt2   (1:nmax,:),&
                  green2 (1:nmax,:),&
		  chil2  (1:nmax,:),&
		  rstpar2(1:nmax,:,:),&
		  topt2  (1:nmax,:),&
		  tll2   (1:nmax,:),&
		  tu2    (1:nmax,:),&
                  defac2 (1:nmax,:),&
		  ph12   (1:nmax,:),&
		  ph22   (1:nmax,:))

      DO i=1,nmax
        ntyp=itype(i)
	topt2	(i,2) =   topt  (ntyp,2)
	tll2	(i,2) =   tll   (ntyp,2)
	tu2	(i,2) =   tu    (ntyp,2)
	defac2  (i,2) =   defac (ntyp,2)
	ph12	(i,2) =   ph1   (ntyp,2)
	ph22	(i,2) =   ph2   (ntyp,2)
	DO  nn = 1, 3
          rstpar2 (i,2,nn)= rstpar(ntyp,2,nn)
        END DO
	zlt2   (i,2) =  zlt   (ntyp,month2(i),2)
	green2 (i,2) =  green (ntyp,month2(i),2)
	chil2  (i,2) =  chil  (ntyp,2)
        vcover (i,2) =  xcover(ntyp,month2(i),2)
        rbc    (i)   =  xbc   (ntyp,month2(i))
        rdc    (i)   =  xdc   (ntyp,month2(i))
        z0x    (i)   =  x0x   (ntyp,month2(i))
        d      (i)   =  xd    (ntyp,month2(i))
        td     (i)   =  tdm   (i)
        tg     (i)   =  tgm   (i)
        tc     (i)   =  tcm   (i)
        capac  (i,1) =  capacm(i,1)
        capac  (i,2) =  capacm(i,2)
        w      (i,1) =  wm    (i,1)
        w      (i,2) =  wm    (i,2)
        w      (i,3) =  wm    (i,3)
      END DO
      ncount=0
      DO i=1,ncols
        IF(imask(i).GE.1) THEN
           ncount=ncount+1
           psur(i)=gps(i)
           tm  (i)=gt (i,1)
           qm  (i)=gq (i,1)
           um  (i)=gu (i,1)/ sin( colrad(i))
           vm  (i)=gv (i,1)/ sin( colrad(i))
           cosz(i)=zenith(i)
        END IF
      END DO
      cpsy=cp/(hl*epsfac)
      rbyg=gasr/grav*delsig(1)*0.5e0
      r100=100.0e0 /gasr
      DO i=1,nmax
        psy   (i)=cpsy*psur(i)
        bps   (i)=sigki(1)
        psb   (i)=psur(i)*delsig(1)
        em    (i)=qm(i)*psur(i)/(epsfac+qm(i))
        dzm   (i)=rbyg*tm(i)
        rhoair(i)=r100*psur(i)/tm(i)
        rcp   (i)=cp*rhoair(i)
      END DO
      ncount=0
      DO i=1,ncols
        IF(imask(i).GE.1) THEN
           ncount=ncount+1
           gmt(i,1)=tmtx(i,1,1)
           gmt(i,2)=tmtx(i,1,2)
           gmt(i,3)=tmtx(i,1,3)
           gmq(i,1)=qmtx(i,1,1)
           gmq(i,2)=qmtx(i,1,2)
           gmq(i,3)=qmtx(i,1,3)
           gmu(i,1)=umtx(i,1,1)
           gmu(i,2)=umtx(i,1,2)
           gmu(i,3)=umtx(i,1,3)
           gmu(i,4)=umtx(i,1,4)
        END IF
      END DO
  END IF
  !
  !     surface albedo (vis/nir and beam/diffuse)
  !     extinction coefficients
  !
  IF(nmax.GE.1) THEN
      CALL radalb ( &
              ncols    , &! INTENT(IN   ) :: ncols
	      month2   , &! INTENT(INOUT) :: mon   (ncols)  	      
	      nmax     , &! INTENT(IN   ) :: nmax		      
	      itype    , &! INTENT(IN   ) :: itype (ncols)	      
	      vcover   , &! INTENT(IN   ) :: vcover(ncols,icg)	      
	      tc       , &! INTENT(IN   ) :: tc    (ncols)	      
	      tg       , &! INTENT(IN   ) :: tg    (ncols)	      
	      capac    , &! INTENT(IN   ) :: capac (ncols,2)	      
	      satcap   , &! INTENT(OUT  ) :: satcap(ncols,icg)	      
	      extk     , &! INTENT(OUT  ) :: extk  (ncols,icg,iwv,ibd) 
	      radfac   , &! INTENT(OUT  ) :: radfac(ncols,icg,iwv,ibd) 
	      closs    , &! INTENT(OUT  ) :: closs (ncols)	      
	      gloss    , &! INTENT(OUT  ) :: gloss (ncols)	      
	      thermk   , &! INTENT(OUT  ) :: thermk(ncols)	      
	      p1f      , &! INTENT(OUT  ) :: p1f   (ncols)	      
	      p2f      , &! INTENT(OUT  ) :: p2f   (ncols)	      
	      zlwup    , &! INTENT(OUT  ) :: zlwup (ncols)	      
	      salb     , &! INTENT(OUT  ) :: salb  (ncols,2,2)	      
	      tgeff    , &! INTENT(OUT  ) :: tgeff (ncols)	      
	      cosz     , &! INTENT(IN   ) :: cosz  (ncols)	      
	      nsx      , &! INTENT(IN   ) :: nsx 		      
	      zlt2       )! INTENT(in   ) :: zlt2  (ncols,icg)
  END IF
  ncount=0
  DO i=1,ncols
      IF(imask(i).GE.1) THEN
        ncount=ncount+1
        avisb(i)=salb(i,1,1)
        avisd(i)=salb(i,1,2)
        anirb(i)=salb(i,2,1)
        anird(i)=salb(i,2,2)
      ELSE IF(ABS(tsea(i)).GE.271.16e0 +0.01e0) THEN
        f=MAX(zenith(i),0.0e0 )
        ocealb=0.12347e0 +f*(0.34667e0+f*(-1.7485e0 + &
             f*(2.04630e0 -0.74839e0 *f)))
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
    !  tsurf(i)=ABS(tsea(i))
    !  qsurf(i)=0.622e0*EXP(21.65605e0 -5418.0e0 /tsurf(i))/gps(i)
  END DO
!  DO i=1,ncols
!      tsurf(i)=ABS(tsea(i))
!      qsurf(i)=0.622e0*EXP(21.65605e0 -5418.0e0 /tsurf(i))/gps(i)
!  END DO
  !
  !     sw/ir heating/cooling rate
  !     downward shortwave radiation at the surface
  !     downward  longwave radiation at the surface
  !
  jlat  = latco




  tice01=tice+0.01
  DO i=1,ncols
      tsurf(i)=ABS(tsea(i))
      IF(initlz.GE.0.AND.kt.EQ.0.AND.ktm.EQ.-1) THEN
        IF(tsea(i).GT.0.0e0.OR. &
             (tsea(i).LT.0.0e0.AND.ABS(tsea(i)).LT.tice01) ) THEN
           tsurf(i)=gt(i,1)
        END IF
        IF(tsea(i).GT.0.0e0.AND.sheleg(i).GT.0.0e0) THEN
           tsurf(i)=MIN(tf  ,tsurf(i))
        END IF
        IF(tsea(i).LT.0.0e0.AND.ABS(tsea(i)).LT.tice01) THEN
           tsurf(i)=MIN(tice,tsurf(i))
        END IF
      END IF
      qsurf(i)=0.622e0*EXP(21.65605e0 -5418.0e0 /tsurf(i))/gps(i)
  END DO
  !
  !     inalb=0:albedo is calculated outside the subr.spmrad.
  !     inalb=2:albedo is calculated outside the subr.spmrad.
  !
  inalb=2
  IF(irad.NE.'YES ')go to 608
  CALL rqvirt(tenth*gps , &! INTENT(in   ) :: pstar (ncols)
              gqq       , &! INTENT(in   ) :: qin   (ncols,kmax) 
	      gtmp      , &! INTENT(inout) :: t     (ncols,kmax)
	      q         , &! INTENT(out  ) :: qout  (ncols,kmax)
	      relhum    , &! INTENT(out  ) :: relhum(ncols,kmax)
	      sig       , &! INTENT(in   ) :: sl    (kmax)
	      ncols     , &! INTENT(IN   ) :: ncols
              kmax        )! INTENT(IN   ) :: kmax
  !
  !     virtual temperature correction for temperature used in radiation
  !     and setting minimum specific humidity
  !
  !     pcb     surface pressure in cb
  !     tsurf   earth's surface temperature used for radiation
  !     for the first time step when ground temperature is not yet
  !     computed (this is done by subr.tsinit ),
  !     gt(k=1)(with some correction below) is used.
  !
  icld=1
  IF(iccon.EQ.'ARA ')icld=3

  CALL spmrad (&
       idatec,&! INTENT(in   ) :: id(4)
       latco ,&! INTENT(in   ) :: jlat 
       colrad,&! INTENT(in   ) :: colrad(ncols)
       gps   ,&! INTENT(in   ) :: gps	(ncols)  
       gtmp  ,&! INTENT(in   ) :: gtmp  (ncols,kmax)
       q     ,&! INTENT(in   ) :: gwv	(ncols,kmax)   
       tsurf ,&! INTENT(in   ) :: gtg	(ncols)   
       relhum,&! INTENT(in   ) :: grh	(ncols,kmax) 
       omg   ,&! INTENT(in   ) :: omg	(ncols,kmax) 
       clr   ,&! INTENT(out  ) :: clr	(ncols,kmax) 
       htr   ,&! INTENT(out  ) :: htr	(ncols,kmax)
       avisd ,&! INTENT(in   ) :: avisd (ncols)
       anird ,&! INTENT(in   ) :: anird (ncols)
       avisb ,&! INTENT(inout) :: avisb (ncols)
       anirb ,&! INTENT(inout) :: anirb (ncols)
       dlwbot,&! INTENT(out  ) :: dlwbot(ncols)
       ulwtop,&! INTENT(out  ) :: ulwtop(ncols)
       dswtop,&! INTENT(out  ) :: dswtop(ncols)
       sig   ,&! INTENT(in   ) :: sig   (kmax)
       sigml ,&! INTENT(in   ) :: sigml (kmax)
       rvisd ,&! INTENT(out  ) :: rvisd (ncols)
       rnird ,&! INTENT(out  ) :: rnird (ncols)   
       rvisb ,&! INTENT(out  ) :: rvisb (ncols)
       rnirb ,&! INTENT(out  ) :: rnirb (ncols)
       inalb ,&! INTENT(in   ) :: inalb 
       trint ,&! INTENT(in   ) :: trint 
       swint ,&! INTENT(in   ) :: swint 
       ultclr,&! INTENT(out  ) :: ultclr(ncols)
       rsclr ,&! INTENT(out  ) :: rsclr (ncols)
       dlwclr,&! INTENT(out  ) :: dlwclr(ncols)
       htrc  ,&! INTENT(out  ) :: htrc  (ncols,kmax)
       rvisdc,&! INTENT(out  ) :: rvisdc(ncols)
       rnirdc,&! INTENT(out  ) :: rnirdc(ncols)
       rvisbc,&! INTENT(out  ) :: rvisbc(ncols)
       rnirbc,&! INTENT(out  ) :: rnirbc(ncols)
       yvisdc,&! INTENT(out  ) :: yvisdc(ncols)
       ynirdc,&! INTENT(out  ) :: ynirdc(ncols)
       yvisbc,&! INTENT(out  ) :: yvisbc(ncols)
       ynirbc,&! INTENT(out  ) :: ynirbc(ncols)
       icld  ,&! INTENT(in   ) :: icld
       zenith,&! INTENT(out  ) :: cosz  (ncols)
       yvisd ,&! INTENT(out  ) :: yvisd (ncols)
       ynird ,&! INTENT(out  ) :: ynird (ncols)
       yvisb ,&! INTENT(out  ) :: yvisb (ncols)
       ynirb ,&! INTENT(out  ) :: ynirb (ncols)
       yswtop,&! INTENT(out  ) :: yswtop(ncols)
       xvisdc,&! INTENT(out  ) :: xvisdc(ncols)
       xnirdc,&! INTENT(out  ) :: xnirdc(ncols)
       xvisbc,&! INTENT(out  ) :: xvisbc(ncols)
       xnirbc,&! INTENT(out  ) :: xnirbc(ncols)
       xvisd ,&! INTENT(out  ) :: xvisd (ncols)
       xnird ,&! INTENT(out  ) :: xnird (ncols)
       xvisb ,&! INTENT(out  ) :: xvisb (ncols)
       xnirb ,&! INTENT(out  ) :: xnirb (ncols)
       cldsav,&! INTENT(out  ) :: cld3  (ncols)
       rs    ,&! INTENT(out  ) :: rs    (ncols)
       ifday ,&! INTENT(in   ) :: ifday
       tod   ,&! INTENT(in   ) :: tod
       idate ,&! INTENT(in   ) :: idate(4)
       cldtot,&! INTENT(out  ) :: cldtd (ncols,kmax)
       cldinv,&! INTENT(out  ) :: cldin (ncols,kmax)
       cldsat,&! INTENT(out  ) :: cstcin(ncols,kmax)
       cldcon,&! INTENT(out  ) :: cldcon(ncols,kmax)
       cldson,&! INTENT(out  ) :: cldson(ncols,kmax)
       clwd  ,&! INTENT(out  ) :: clwd  (ncols,kmax)
       emisd ,&! INTENT(out  ) :: emisd (ncols,kmax)
       taud  ,&! INTENT(out  ) :: taud  (ncols,kmax)
       tsea  ,&! INTENT(in   ) :: tsea  (ncols)
       yrl   ,&! INTENT(in   ) :: yrl
       monl  ,&! INTENT(in   ) :: monl(12)
       convts,&! INTENT(in   ) :: convts(ncols)
       convcs,&! INTENT(in   ) :: convcs(ncols)
       convbs,&! INTENT(in   ) :: convbs(ncols)
       mxrdcc,&! intent(in   ) :: mxrdcc
       lcnvl ,&! intent(in   ) :: lcnvl  
       lthncl,&! intent(in   ) :: lthncl 
       convc ,&! intent(in   ) :: convc (ncols)
       convt ,&! intent(in   ) :: convt (ncols)
       convb ,&! intent(in   ) :: convb (ncols)
       first ,&! INTENT(in   ) :: first
       nfprt ,&! INTENT(in   ) :: nfprt     
       nferr ,&! INTENT(in   ) :: nferr     
       ifprt ,&! INTENT(in   ) :: ifprt(100)
       dodia ,&! INTENT(in   ) :: dodia(ndavl) 
       ndavl ,&! INTENT(IN   ) :: ndavl 
       co2val,&! INTENT(in   ) :: co2val
       delt  ,&! INTENT(in   ) :: delt
       nfin0 ,&! INTENT(in   ) :: nfin0
       nfin1 ,&! INTENT(in   ) :: nfin1
       nfcnv0,&! INTENT(in   ) :: nfcnv0 
       nls   ,&! INTENT(IN   ) :: nls  
       nlcs  ,&! INTENT(IN   ) :: nlcs 
       ncols ,&! INTENT(IN   ) :: ncols 
       jmax  ,&! INTENT(IN   ) :: jmax  
       kmax  ,&! INTENT(IN   ) :: kmax  
       lonrad,&! INTENT(in   ) :: lonrad(ncols)
       cos2d ,&! INTENT(IN   ) :: cos2d (ncols)
       intcosz)! INTENT(IN   ) :: intcosz
  !
  !     save cloud radiation fields from initial spmrad call for 
  !     future inclusion into diagnostics accumulators
  !
  IF(cloudr.AND.jdt.EQ.1.AND.ktm.EQ.-1.AND.kt.EQ.0)THEN
      WRITE(nfcldr,rec=latco)cldtot,cldinv,cldsat,cldcon,cldson,clwd, &
          emisd,taud
  END IF
  !
  !     retrieve cloud radiation fields from initial spmrad call for
  !     inclusion into diagnostics accumulators
  !
  IF(cloudr.AND.jdt.EQ.2.AND.ktm.NE.-1.AND.kt.EQ.0)THEN
      READ(nfcldr,rec=latco)cldtot,cldinv,cldsat,cldcon,cldson,clwd, &
          emisd,taud
  END IF
  IF(ghl_local) THEN
     IF (dogrh(ngctd,latco)) CALL StoreGridHistory (cldtot, ngctd, latco)
  END IF
		 
  IF(cloudr.AND.(ABS( MOD((tod-delt)/f3600+0.03125e0,swint)) &
       .LE.0.0625e0.OR.(jdt.EQ.2.AND.ktm.NE.-1.AND.kt.EQ.0)))THEN
      !
      !     radint is the number of time steps between short wave radiation
      !     calculation, hence the cloud radiation interval
      !
      radint=swint*f3600/delt
      !
      !     use relhum as scratch pad
      !
      IF(dodia(ndctd))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*cldtot(i,k)
           END DO
        END DO
        CALL updia(relhum,ndctd,latco)
      END IF
      IF(dodia(ndcin))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*cldinv(i,k)
           END DO
        END DO
        CALL updia(relhum,ndcin,latco)
      END IF
      IF(dodia(ndcst))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*cldsat(i,k)
           END DO
        END DO
        CALL updia(relhum,ndcst,latco)
      END IF             

      IF(dodia(ndccn))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*cldcon(i,k)
           END DO
        END DO
        CALL updia(relhum,ndccn,latco)
      END IF
      IF(dodia(ndcsn))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*cldson(i,k)
           END DO
        END DO
        CALL updia(relhum,ndcsn,latco)
      END IF
      IF(dodia(ndclw))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*clwd(i,k)
           END DO
        END DO
        CALL updia(relhum,ndclw,latco)
      END IF
      IF(dodia(ndsod))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*taud(i,k)
           END DO
        END DO
        CALL updia(relhum,ndsod,latco)
      END IF
  END IF
  IF(cloudr.AND.(ABS( MOD((tod-delt)/f3600+0.03125e0,trint)) &
       .LE.0.0625e0.OR.(jdt.EQ.2.AND.ktm.NE.-1.AND.kt.EQ.0)))THEN
      !
      !     radint is the number of time steps between long wave radiation
      !     calculation, hence the cloud radiation interval
      !
      radint=trint*f3600/delt
      !
      !     use relhum as scratch pad
      !
      IF(dodia(ndemi))THEN
        DO k=1,kmax
           DO i=1,ncols
              relhum(i,k)=radint*emisd(i,k)
           END DO
        END DO
        CALL updia(relhum,ndemi,latco)
      END IF
  END IF
  !
  !     td  initialization
  !     there are two problems of this scheme
  !     over antactica tg seems to be too high compared with control tg0
  !     (or contol tg0 may be too low)
  !     in   siberia   td seems to be too low
  !     modification is recommended
  !
  IF(initlz.GE.0.AND.ktm.EQ.-1.AND.kt.EQ.0.AND.nmax.GE.001)THEN
      DO ind=1,2
        ncount=0
        DO i=1,ncols
           IF(imask(i).GE.1) THEN
              ncount=ncount+1
              IF(ind.EQ.1) THEN
                 !
                 !     night
                 !
                 radn(i,1,1)=0.0e0
                 radn(i,1,2)=0.0e0
                 radn(i,2,1)=0.0e0
                 radn(i,2,2)=0.0e0
                 cosz(i)    =0.0e0
              ELSE
                 !
                 !     noon
                 !
                 radn(i,1,1)=xvisb (i)
                 radn(i,1,2)=xvisd (i)
                 radn(i,2,1)=xnirb (i)
                 radn(i,2,2)=xnird (i)
                 cosz(i)    =cos2(i)
              END IF
              radn(i,3,1)=0.0e0
              radn(i,3,2)=dlwbot(i)
              !
              !     precipitation
              !
              ppl (i)    =0.0e0
              ppc (i)    =0.0e0
           END IF
        END DO
        DO itr=1,5
           CALL radalb( &
	       ncols ,month2   ,nmax  ,itype ,vcover,tc    , &
	       tg    ,capac ,satcap,extk  ,radfac,closs ,gloss ,thermk,p1f   , &
	       p2f   ,zlwup ,salb  ,tgeff ,cosz  ,nsx  , zlt2)
	
	
        CALL fysiks(vcover, z0x  , d    , rdc  , rbc  , z0   , ndt  , &
		    latitu, nfprt, bps  ,psb   ,dzm   ,em    ,gmt   ,gmq   , &
		    gmu   ,cu    ,cuni  ,ctni  ,ustar ,cosz  ,sinclt2,rhoair, &
		    psy   ,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    , ect  , &
		    eci   , egt  , egi  , egs  , ec   , eg   , hc   , hg   , &
		    ecidif,egidif,ecmass,egmass,etmass, hflux, chf  , shf  , &
		    fluxef, roff , drag ,ra    , rb   , rd   , rc   , rg   , &
		    tcta  , tgta , ta   , ea   , etc  , etg  , btc  , btg  , &
		    u2    , radt , par  , pd   , rst  ,rsoil ,phroot, hrr  , &
		    phsoil, cc   , cg   ,satcap, snow , dtc  , dtg  , dtm  , &
		    dqm   , stm  , extk ,radfac, closs,gloss ,thermk, p1f  , &
		    p2f   , tc   , tg   , td   , capac, w    ,  qm  , tm   , &
		    um    , vm   , psur , ppc  , ppl  , radn ,itype , dtc3x, &
		    month2, nmax , ncols,zlt2  ,green2,chil2 ,rstpar2,topt2, &
		    tll2  ,tu2   , defac2,ph12  ,ph22  )
	   
	   
           ncount=0
           DO i=1,ncols
              IF(imask(i).GE.1) THEN
                 ncount=ncount+1
                 tm (i  )=gt  (i,1)
                 qm (i  )=gq  (i,1)
                 gmt(i,1)=tmtx(i,1,1)
                 gmt(i,2)=tmtx(i,1,2)
                 gmt(i,3)=tmtx(i,1,3)
                 gmq(i,1)=qmtx(i,1,1)
                 gmq(i,2)=qmtx(i,1,2)
                 gmq(i,3)=qmtx(i,1,3)
                 gmu(i,1)=umtx(i,1,1)
                 gmu(i,2)=umtx(i,1,2)
                 gmu(i,3)=umtx(i,1,3)
                 gmu(i,4)=umtx(i,1,4)
              END IF
           END DO
        END DO
        DO i=1,nmax
           capac(i,1)=capacm(i,1)
           capac(i,2)=capacm(i,2)
           w    (i,1)=wm(i,1)
           w    (i,2)=wm(i,2)
           w    (i,3)=wm(i,3)
           td   (i)  =tdm(i)
           tc   (i)  =tcm(i)
           IF(ind.EQ.1) THEN
              tmin (i) =tg (i)
           ELSE
              tmax (i) =tg (i)
           END IF
           tg   (i) =tgm(i)
        END DO
      END DO
      DO i=1,nmax
        td   (i) =0.9*0.5*(tmax(i)+tmin(i))+0.1*tdm(i)
        tdm  (i) =td(i)
        td0  (i) =td(i)
      END DO
      !
      !     this is a start of equilibrium tg,tc comp.
      !     
      ncount=0
      DO i=1,ncols
        IF(imask(i).GE.1) THEN
           ncount=ncount+1
           cosz(i)    =zenith(i)
        END IF
      END DO
      DO i=1,nmax
        IF(cosz(i).LT.0.0e0) THEN
           tgm  (i)  =tmin(i)
           tg0  (i)  =tmin(i)
        END IF
      END DO
      CALL radalb ( &
              ncols ,month2,nmax  ,itype ,vcover,tc    , &
	      tg    ,capac ,satcap,extk  ,radfac,closs ,gloss ,thermk,p1f   , &
	      p2f   ,zlwup ,salb  ,tgeff ,cosz  ,nsx   ,zlt2)
  END IF
608 CONTINUE
  !
  !     this is for radiation interpolation
  !
  intstp=3600.0*swint/delt+0.01
  modstp=MOD(jdt-1,intstp)
  fstp  =float(modstp)/float(intstp)
  IF(jdt.LE.2.AND.kt.EQ.0) fstp=0.0
  fstp1 =1.0-fstp
  DO i=1,ncols
      xvisd (i)=fstp1*rvisd (i)+fstp*yvisd (i)
      xvisb (i)=fstp1*rvisb (i)+fstp*yvisb (i)
      xnird (i)=fstp1*rnird (i)+fstp*ynird (i)
      xnirb (i)=fstp1*rnirb (i)+fstp*ynirb (i)
      xvisdc(i)=fstp1*rvisdc(i)+fstp*yvisdc(i)
      xvisbc(i)=fstp1*rvisbc(i)+fstp*yvisbc(i)
      xnirdc(i)=fstp1*rnirdc(i)+fstp*ynirdc(i)
      xnirbc(i)=fstp1*rnirbc(i)+fstp*ynirbc(i)
      xswtop(i)=fstp1*dswtop(i)+fstp*yswtop(i)
  END DO
  IF(modstp.EQ.intstp-1) THEN
      DO i=1,ncols
        rvisd (i)=yvisd (i)
        rvisb (i)=yvisb (i)
        rnird (i)=ynird (i)
        rnirb (i)=ynirb (i)
        rvisdc(i)=yvisdc(i)
        rvisbc(i)=yvisbc(i)
        rnirdc(i)=ynirdc(i)
        rnirbc(i)=ynirbc(i)
        dswtop(i)=yswtop(i)
      END DO
  END IF
  IF(nmax.GE.1) THEN
      ncount=0
      DO i=1,ncols
        IF(imask(i).GE.1) THEN
           ncount=ncount+1
           !
           !     this is for radiation interpolation
           !
           IF(cosz(i).GE.0.01746e0 ) THEN
              radn(i,1,1)=xvisb (i)
              radn(i,1,2)=xvisd (i)
              radn(i,2,1)=xnirb (i)
              radn(i,2,2)=xnird (i)
           ELSE
              radn(i,1,1)=0.0e0
              radn(i,1,2)=0.0e0
              radn(i,2,1)=0.0e0
              radn(i,2,2)=0.0e0
           END IF
           radn(i,3,1)=0.0e0
           radn(i,3,2)=dlwbot(i)
           !
           !     precipitation
           !
           ppl (i)    =ppli  (i)
           ppc (i)    =ppci  (i)
        END IF
      END DO
  END IF
  !
  !     surface parameterization
  !     simple biosphere parameterization
  !
  IF(nmax.GE.1) THEN
  
      CALL fysiks(vcover, z0x  , d    , rdc  , rbc  , z0   ,ndt  , &
		 latitu, nfprt, bps  ,psb   ,dzm   ,em    ,gmt   ,gmq   , &
		 gmu   ,cu    ,cuni  ,ctni  ,ustar ,cosz  ,sinclt2,rhoair, &
		 psy   ,rcp   ,wc    ,wg    ,fc    ,fg    ,hr    , ect  , &
		 eci   , egt  , egi  , egs  , ec   , eg   , hc   , hg   , &
		 ecidif,egidif,ecmass,egmass,etmass,hflux , chf  , shf  , &
		 fluxef, roff , drag ,ra    , rb   , rd   , rc   , rg   , &
		 tcta  , tgta , ta   , ea   , etc  , etg  , btc  , btg  , &
		 u2    , radt , par  , pd   , rst  ,rsoil ,phroot, hrr  , &
		 phsoil, cc   , cg   ,satcap, snow , dtc  , dtg  , dtm  , &
		 dqm   , stm  , extk ,radfac, closs,gloss ,thermk, p1f  , &
		 p2f   , tc   , tg   , td   , capac, w    , qm   , tm   , &
		 um    , vm   , psur , ppc  , ppl  , radn ,itype , dtc3x, &
		 month2, nmax , ncols,zlt2  ,green2,chil2 ,rstpar2,topt2, &
		 tll2  ,tu2   , defac2,ph12  ,ph22)
     
      ncount=0
      DO i=1,ncols
        IF(imask(i).GE.1) THEN
           ncount=ncount+1
           tmtx(i,1,3)=gmt(i,3)
           qmtx(i,1,3)=gmq(i,3)
           umtx(i,1,3)=gmu(i,3)
           umtx(i,1,4)=gmu(i,4)
           tsea(i)    =tgeff(i)
        END IF
      END DO
  END IF
  !
  !     sea or sea ice
  !
  mmax=ncols-nmax+1
  IF(mmax.GE.1) THEN
      DO i=1,ncols
        !
        !     this is for radiation interpolation
        !
        slrad(i)=dlwbot(i)+(1.0-avisd(i))*xvisd(i) &
             +(1.0-avisb(i))*xvisb(i) &
             +(1.0-anird(i))*xnird(i) &
             +(1.0-anirb(i))*xnirb(i)
        slrad(i)=-cnvfac*slrad(i)
        xsea (i)=        tseam(i)
      END DO
      
      CALL seasfc( &
              tmtx  ,umtx  ,qmtx  ,kmax  ,kmax  ,slrad ,tsurf ,qsurf , &
              gu    ,gv    ,gt    ,gq    ,gps   ,xsea  ,dtc3x ,sin(colrad), &
	      sigki ,delsig,sens  ,evap  ,umom  ,vmom  ,rmi   ,rhi   , &
	      cond  ,stor  ,zorl  ,rnet  ,ncols)
      
      DO i=1,ncols
        IF(tsea(i).LE.0.0e0.AND.tsurf(i).LT.tice+0.01e0 ) THEN
           IF(intg.EQ.2) THEN
              IF(istrt.EQ.0) THEN
                 tseam(i)=filta*tsea(i)+epsflt*(tseam(i)+xsea(i))
              END IF
              tsea (i)=      xsea(i)
           ELSE
              tsea (i)=xsea(i)
              tseam(i)=xsea(i)
           END IF
        END IF
        IF(tsea(i).LT.0.0e0.AND.tsurf(i).GE.tice+0.01e0) THEN
           tseam(i) = tsea(i)
        ENDIF
      END DO
  END IF
  !
  !     yamada-mellor pbl parameterization
  !     ( solving the matrices from bottom to top )
  !     
  CALL ympbl1 (tmtx,qmtx,umtx,ncols,kmax,kmax)

  DO k=1,kmax
      DO i=1,ncols
        !gyu(i,k)=gyu(i,k)-umtx(i,k,4)
        !gyv(i,k)=gyv(i,k)+umtx(i,k,3)     
	gyu(i,k)  = gyu(i,k) + umtx(i,k,3)
        gyv(i,k)  = gyv(i,k) + umtx(i,k,4)
      END DO
  END DO
 
  DO k=1,kmax
      DO i=1,ncols
        gqd(i,k)=gqd(i,k)+qmtx(i,k,3)
      END DO
  END DO
  !     
  !     tmtx(i,k,3) virtual temperature tendency due to vertical diffusion
  !
  DO k=1,kmax
      DO i=1,ncols
        tmtx(i,k,3) = (one + delq * gq(i,k)) * tmtx(i,k,3) + &
                       delq * gt(i,k) * qmtx(i,k,3)
        gtd(i,k)=gtd(i,k) + tmtx(i,k,3)
        gtd(i,k)=gtd(i,k) + (one+delq*gq(i,k)) * (htr(i,k)+clr(i,k))
      END DO
  END DO
  IF(kmax.LT.kmax) THEN
      DO k=kmax+1,kmax
        DO i=1,ncols
           gtd(i,k)=gtd(i,k) + htr(i,k)+clr(i,k)
           gtd(i,k)=gtd(i,k) + tmtx(i,k,3)
        END DO
      END DO
  END IF
  !
  !     pointwise diagnostics
  !
  ncount=0
  DO i=1,ncols
      IF(imask(i).GE.1) THEN
        ncount=ncount+1
        fmom= rhoair(i)*cu(i)*ustar(i)
        umom(i)=fmom*um(i)
        vmom(i)=fmom*vm(i)
	!JPB
	strcut=10.0
	umom(i)=MIN(strcut,umom(i))
	vmom(i)=MIN(strcut,vmom(i))
	umom(i)=MAX(-strcut,umom(i))
	vmom(i)=MAX(-strcut,vmom(i))
	!JPB
      END IF
  END DO
  IF (ghl_local) THEN
      IF (dogrh(nglwh,latco) ) CALL StoreGridHistory (clr,    nglwh,  latco)
      IF (dogrh(ngswh,latco) ) CALL StoreGridHistory (htr,    ngswh,  latco)
      IF (dogrh(ngulwt,latco)) CALL StoreGridHistory (ulwtop, ngulwt, latco)
      IF (dogrh(ngdlwb,latco)) CALL StoreGridHistory (dlwbot, ngdlwb, latco)
      IF (dogrh(ngcsz,latco)) CALL StoreGridHistory (zenith, ngcsz,  latco)
      IF (dogrh(ngdswt,latco)) CALL StoreGridHistory (xswtop, ngdswt, latco)
      IF (dogrh(ngsvb,latco) ) CALL StoreGridHistory (xvisb,  ngsvb,  latco)
      IF (dogrh(ngsvd,latco) ) CALL StoreGridHistory (xvisd,  ngsvd,  latco)
      IF (dogrh(ngsnb,latco) ) CALL StoreGridHistory (xnirb,  ngsnb,  latco)
      IF (dogrh(ngsnd,latco) ) CALL StoreGridHistory (xnird,  ngsnd,  latco)
      IF (dogrh(ngavb,latco) ) CALL StoreGridHistory (avisb,  ngavb,  latco)
      IF (dogrh(ngavd,latco) ) CALL StoreGridHistory (avisd,  ngavd,  latco)
      IF (dogrh(nganb,latco) ) CALL StoreGridHistory (anirb,  nganb,  latco)
      IF (dogrh(ngand,latco) ) CALL StoreGridHistory (anird,  ngand,  latco)
      IF (dogrh(ngcnd,latco) ) CALL StoreGridHistory (cond,   ngcnd,  latco)
      IF (dogrh(ngstr,latco) ) CALL StoreGridHistory (stor,   ngstr,  latco)
  END IF
  !
  !     diagnostics of radiation
  !
  IF( (kt.NE.0) .OR. (jdt.NE.1) ) THEN
      IF(dodia(nddtlw))   &
        CALL updia(clr,nddtlw,latco)
		  
      IF(dodia(nddtsw)) &
       CALL updia(htr,nddtsw,latco)
		  
      IF(dodia(ndcld)) &
       CALL  updia(cldsav,ndcld,latco)
		  
      IF (ghl_local) THEN
         IF(dogrh(ngcld,latco)) CALL StoreGridHistory(cldsav, ngcld, latco)
      END IF

       uswtop=0.0   ! CALL reset(uswtop,ncols)
       uswtpc=0.0   ! CALL reset(uswtpc,ncols)
      
      DO k=1,kmax
         rfac=100.0e0*cp*delsig(k)/grav
         DO i=1,ncols
            uswtop(i)=uswtop(i)+rfac*gps(i)*htr(i,k)
            uswtpc(i)=uswtpc(i)+rfac*gps(i)*htrc(i,k)
         END DO
      END DO
      DO i=1,ncols
         uswtop(i)=xswtop(i) &
              -(1.0-avisd(i))*xvisd(i) &
              -(1.0-avisb(i))*xvisb(i) &
              -(1.0-anird(i))*xnird(i) &
              -(1.0-anirb(i))*xnirb(i) &
              -uswtop(i)
         uswtpc(i)=xswtop(i) &
              -(1.0-avisd(i))*xvisdc(i) &
              -(1.0-avisb(i))*xvisbc(i) &
              -(1.0-anird(i))*xnirdc(i) &
              -(1.0-anirb(i))*xnirbc(i) &
              -uswtpc(i)
      END DO
      IF(dodia(ndswea))THEN
         DO i=1,ncols
            bfr1(i)=xswtop(i)-uswtop(i)
         END DO
         CALL updia(bfr1,ndswea,latco)
      END IF
      IF(dodia(ndseac))THEN
         DO i=1,ncols
            bfr1(i)=xswtop(i)-uswtpc(i)
         END DO
         CALL updia(bfr1,ndseac,latco)
		  
      END IF
      IF(dodia(ndlwdb)) &
         CALL updia(dlwbot,ndlwdb,latco)
		  
      IF(dodia(ndldbc)) &
        CALL updia(dlwclr,ndldbc,latco)
		  
      IF(dodia(ndlwub))THEN
         DO i=1,ncols
            bfr1(i)=stefan*tsurf(i)**4
         END DO
         CALL updia(bfr1,ndlwub,latco)
		  
      END IF
      IF(dodia(ndlwut)) &
       CALL updia(ulwtop,ndlwut,latco)
		  
      IF(dodia(ndlutc)) &
       CALL updia(ultclr,ndlutc,latco)
		  
      IF(dodia(ndlwnb)) & 
       CALL updia(rs,ndlwnb,latco)
		  
      IF(dodia(ndlnbc)) &
       CALL updia(rsclr,ndlnbc,latco)
		  
      IF(dodia(ndswdb))THEN
         DO i=1,ncols
            bfr1(i)=xvisd(i)+xvisb(i)+xnird(i)+xnirb(i)
         END DO
         CALL updia(bfr1,ndswdb,latco)
		  
      END IF
      IF(dodia(ndsdbc))THEN
        DO i=1,ncols
           bfr1(i)=xvisdc(i)+xvisbc(i)+xnirdc(i)+xnirbc(i)
        END DO
        CALL updia(bfr1,ndsdbc,latco)
      END IF
      IF(dodia(ndswe))THEN
         DO i=1,ncols
            bfr1(i)=(1.0-avisd(i))*xvisd(i)+(1.0-avisb(i))*xvisb(i) &
                 +(1.0-anird(i))*xnird(i)+(1.0-anirb(i))*xnirb(i) 
         END DO
         CALL updia(bfr1,ndswe,latco)
		  
      END IF
      IF(dodia(ndswec))THEN
         DO i=1,ncols
            bfr1(i)=(1.0-avisd(i))*xvisdc(i)+(1.0-avisb(i))*xvisbc(i) &
                    +(1.0-anird(i))*xnirdc(i)+(1.0-anirb(i))*xnirbc(i)
         END DO
         CALL updia(bfr1,ndswec,latco)
      END IF
      IF(dodia(ndswub))THEN
        DO i=1,ncols
           bfr1(i)=avisd(i)*xvisd(i)+avisb(i)*xvisb(i)+ &
                   anird(i)*xnird(i)+anirb(i)*xnirb(i)
        END DO
        CALL updia(bfr1,ndswub,latco)
      END IF
      IF(dodia(ndsubc))THEN
         DO i=1,ncols
            bfr1(i)=avisd(i)*xvisdc(i)+avisb(i)*xvisbc(i)+ &
                    anird(i)*xnirdc(i)+anirb(i)*xnirbc(i)
         END DO
         CALL updia(bfr1,ndsubc,latco)
      END IF
      IF(dodia(ndswdt)) &
       CALL updia(xswtop,ndswdt,latco)
      IF(dodia(ndswut)) &
       CALL updia(uswtop,ndswut,latco)
      IF(dodia(ndsutc)) &
       CALL updia(uswtpc,ndsutc,latco)
      !
      !     diagnostics of pbl parameterization
      !
      dtc3xi=1.0/dtc3x
      IF(ghl_local)THEN
         IF(dogrh(ngrarh,latco))THEN
            ncount=0
            DO i=1,ncols
               IF(imask(i).GE.1) THEN
                  ncount=ncount+1
                  bfrg(i)=ra(i)
               ELSE
                  bfrg(i)=one/rhi(i)
               END IF
            END DO
           CALL StoreGridHistory (bfrg, ngrarh, latco)
         END IF
         IF(dogrh(ngrm,latco))THEN
             ncount=0
             DO i=1,ncols
                IF(imask(i).GE.1) THEN
                   ncount=ncount+1
                   bfrg(i)=0.0e0
                ELSE
                   bfrg(i)=one/rmi(i)
               END IF
             END DO
            CALL StoreGridHistory (bfrg, ngrm, latco)
         END IF
	
        IF(dogrh(ngdrg,latco) ) CALL StoreMaskedGridHistory(drag,      imask, ngdrg,  latco)
        IF(dogrh(ngrd1,latco) ) CALL StoreMaskedGridHistory(radt(1,1), imask, ngrd1,  latco)
        IF(dogrh(ngrd2,latco) ) CALL StoreMaskedGridHistory(radt(1,2), imask, ngrd2,  latco)
        IF(dogrh(ngrb,latco)  ) CALL StoreMaskedGridHistory(rb,        imask, ngrb,   latco)
        IF(dogrh(ngrd,latco)  ) CALL StoreMaskedGridHistory(rd,        imask, ngrd,   latco)
        IF(dogrh(ngrc,latco)  ) CALL StoreMaskedGridHistory(rc,        imask, ngrc,   latco)
        IF(dogrh(ngrg,latco)  ) CALL StoreMaskedGridHistory(rg,        imask, ngrg,   latco)
        IF(dogrh(ngrs,latco)  ) CALL StoreMaskedGridHistory(rsoil,     imask, ngrs,   latco)
        IF(dogrh(ngea,latco)  ) CALL StoreMaskedGridHistory(ea,        imask, ngea,   latco)
        IF(dogrh(ngta,latco)  ) CALL StoreMaskedGridHistory(ta,        imask, ngta,   latco)
        IF(dogrh(nghc,latco)  ) CALL StoreMaskedGridHistory(hc,        imask, nghc,   latco, dtc3xi)
        IF(dogrh(nghg,latco)  ) CALL StoreMaskedGridHistory(hg,        imask, nghg,   latco, dtc3xi)
        IF(dogrh(ngect,latco) ) CALL StoreMaskedGridHistory(ect,       imask, ngect,  latco, dtc3xi)
        IF(dogrh(ngegt,latco) ) CALL StoreMaskedGridHistory(egt,       imask, ngegt,  latco, dtc3xi)
        IF(dogrh(ngeci,latco) ) CALL StoreMaskedGridHistory(eci,       imask, ngeci,  latco, dtc3xi)
        IF(dogrh(ngegi,latco) ) CALL StoreMaskedGridHistory(egi,       imask, ngegi,  latco, dtc3xi)
        IF(dogrh(nglhbs,latco)) CALL StoreMaskedGridHistory(egs,       imask, nglhbs, latco, dtc3xi)
        IF(dogrh(ngchf,latco) ) CALL StoreMaskedGridHistory(chf,       imask, ngchf,  latco)
        IF(dogrh(ngghf,latco) ) CALL StoreMaskedGridHistory(shf,       imask, ngghf,  latco)
        IF(dogrh(ngrof,latco) ) CALL StoreMaskedGridHistory(roff,      imask, ngrof,  latco, 1000.0*dtc3xi)
      END IF
      ncount=0
      DO i=1,ncols
         IF(ghl_local .AND. dogrh(ngulwb,latco))bfrg(i)=stefan*tsurf(i)**4
         IF(dodia(ndintl))bfr1(i)=0.0
         IF(dodia(ndrnof))bfr2(i)=0.0
         IF(dodia(ndta))swrk(i,1)=tsurf(i)
         IF(dodia(ndea))swrk(i,2)=0.0
         IF(dodia(ndlhbs))swrk(i,3)=0.0
         IF(dodia(ndrb))swrk(i,6)=0.0
         IF(dodia(ndrd))swrk(i,7)=0.0
         IF(dodia(ndrg))swrk(i,8)=0.0
         IF(dodia(ndrs))swrk(i,9)=0.0
         IF(imask(i).GE.1) THEN
            ncount=ncount+1
            IF(ghl_local .AND. dogrh(ngulwb,latco))bfrg(i)=zlwup(i)
              sens(i)=(hc(i)+hg(i))*dtc3xi
              evap(i)=(ec(i)+eg(i))*dtc3xi
           IF(dodia(ndintl))bfr1(i)= &
                bfr1(i)+(eci(i)+egi(i))*dtc3xi
           IF(dodia(ndrnof))bfr2(i)=bfr2(i)+roff(i)*dtc3xi*thousd
           IF(dodia(ndta))swrk(i,1)=ta(i)
           IF(dodia(ndea))swrk(i,2)=ea(i)
           IF(dodia(ndlhbs))swrk(i,3)=egs(i)*dtc3xi
           IF(dodia(ndrm))swrk(i,4)=0.0
           IF(dodia(ndrarh))swrk(i,5)=ra(i)
           IF(dodia(ndrb))swrk(i,6)=rb(i)
           IF(dodia(ndrd))swrk(i,7)=rd(i)
           IF(dodia(ndrg))swrk(i,8)=rg(i)
           IF(dodia(ndrs))swrk(i,9)=rs(i)
        ELSE
           IF(ghl_local .AND. dogrh(ngulwb,latco))bfrg(i)=stefan*tsurf(i)**4
           IF(dodia(ndrm))swrk(i,4)=one/rmi(i)
           IF(dodia(ndrarh))swrk(i,5)=one/rhi(i)
        ENDIF
     END DO
     IF(dodia(ndshf)) &
      CALL updia(sens,ndshf,latco)
		  
     IF(dodia(ndlhf)) &
      CALL updia(evap,ndlhf,latco)
		  
     IF(dodia(ndintl)) &
      CALL updia(bfr1,ndintl,latco)
		  
     IF(dodia(ndrnof)) &
      CALL updia(bfr2,ndrnof,latco)
     IF(dodia(ndta)) &
      CALL updia(swrk(1,1),ndta,latco)
     IF(dodia(ndea)) &
      CALL updia(swrk(1,2),ndea,latco)
		  
     IF(dodia(ndlhbs)) &
      CALL updia(swrk(1,3),ndlhbs,latco)
		  
     IF(dodia(ndrm)) &
      CALL updia(swrk(1,4),ndrm,latco)
		  
     IF(dodia(ndrarh)) &
      CALL updia(swrk(1,5),ndrarh,latco)
     IF(dodia(ndrb)) &
      CALL updia(swrk(1,6),ndrb,latco)
     IF(dodia(ndrd)) & 
      CALL updia(swrk(1,7),ndrd,latco)
     IF(dodia(ndrg)) &
      CALL updia(swrk(1,8),ndrg,latco)
     IF(dodia(ndrs)) &
      CALL updia(swrk(1,9),ndrs,latco)
		  
     IF (cdhl(jdt)) THEN
        DO i=1,ncols
           ustr(i)=-umom(i)
           vstr(i)=-vmom(i)
        ENDDO
     ENDIF

     IF(dodia(ndustr))THEN
        DO i=1,ncols
           bfr1(i)=-umom(i)
        END DO
        CALL updia(bfr1,ndustr,latco)
     END IF
     IF(dodia(ndvstr))THEN
        DO i=1,ncols
           bfr2(i)=-vmom(i)
        END DO
        CALL updia(bfr2,ndvstr,latco)
     END IF
     IF(dodia(nddtvd)) &
       CALL updia(tmtx(1,1,3),nddtvd,latco)
		  
     IF(dodia(ndduvd)) &
       CALL updia(umtx(1,1,3),ndduvd,latco)
		  
     IF(dodia(nddvvd)) &
       CALL updia(umtx(1,1,4),nddvvd,latco)
		  
     IF(dodia(nddqvd)) &
      CALL updia(qmtx(1,1,3),nddqvd,latco)
     IF(ghl_local)THEN
        IF(dogrh(ngdtvd,latco)) CALL StoreGridHistory(tmtx(:,:,3), ngdtvd, latco)
        IF(dogrh(ngduvd,latco)) CALL StoreGridHistory(umtx(:,:,3), ngduvd, latco)
        IF(dogrh(ngdvvd,latco)) CALL StoreGridHistory(umtx(:,:,4), ngdvvd, latco)
        IF(dogrh(ngdqvd,latco)) CALL StoreGridHistory(qmtx(:,:,3), ngdqvd, latco)
        IF(dogrh(nguswt,latco)) CALL StoreGridHistory(uswtop,      nguswt, latco)
        IF(dogrh(ngulwb,latco)) CALL StoreGridHistory(bfrg,        ngulwb, latco)
        IF(dogrh(ngustr,latco)) CALL StoreGridHistory(umom,        ngustr, latco)
        IF(dogrh(ngvstr,latco)) CALL StoreGridHistory(vmom,        ngvstr, latco)
        IF(dogrh(ngshf,latco) ) CALL StoreGridHistory(sens,        ngshf,  latco)
        IF(dogrh(nglhf,latco) ) CALL StoreGridHistory(evap,        nglhf,  latco)
     END IF
     IF(dodia(ndtsb)) &
       CALL updia(tsurf,ndtsb,latco)
		  
  END IF
  !
  !     sib time integaration and time filter
  !
  IF(nmax.GE.1) THEN
      CALL sextrp ( &
               td    ,tg    ,tc    ,w     ,capac ,td0   ,tg0   ,tc0   ,w0    , &
	       capac0,tdm   ,tgm   ,tcm   ,wm    ,capacm,istrt ,ncols ,nmax  , &
	       epsflt,intg  ,latco ,imask ,nfprt)
      !
      !     fix soil moisture at selected locations
      !
      DO i=1,nmax
        IF(ssib(i).GT.0.0)THEN
           w0(i,1)=ssib(i)
           w0(i,2)=ssib(i)
           w0(i,3)=ssib(i) 
           wm(i,1)=ssib(i)
           wm(i,2)=ssib(i)
           wm(i,3)=ssib(i)
        END IF
      END DO
  END IF  
  DO k=1,kmax
      DO i=1,ncols
        gt(i,k)=gt(i,k)*(1.0e0+0.608e0*gq(i,k))
      END DO
  END DO
 END SUBROUTINE physcs
!***************************************************************************

REAL FUNCTION julday (id,tod,monl)
 IMPLICIT NONE
  INTEGER, INTENT(IN   ) :: id(4)
  INTEGER, INTENT(IN   ) :: monl(12)
  REAL   , INTENT(IN   ) :: tod
  REAL                   :: xday
  REAL                   :: f3600 = 3.6e3  
  REAL                   :: yrl   =365.25e0  
  INTEGER                :: monday(12)
  INTEGER                :: m
!
! compute the julian day from a normal date
!
    monday(1)=0
    DO m=2,12
       monday(m)=monday(m-1)+monl(m-1)
    END DO
   
    xday=id(1)*f3600
    xday=xday+MOD(tod,f3600)
    xday=monday(id(2))+id(3)+xday/86400.0
    xday=xday-MOD(id(4)+3,4)*.25
    IF(MOD(id(4),4).EQ.0.AND.id(2).GT.2)xday=xday+1.0e0
    julday= MOD(xday-1.0e0,yrl)
    
END FUNCTION julday
!*****************************************************************************
