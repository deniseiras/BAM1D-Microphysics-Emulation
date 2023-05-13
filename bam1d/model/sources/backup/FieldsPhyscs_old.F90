MODULE FieldsPhyscs
  USE Constants, ONLY :	  &
       r8,i8,r4,i4
      
  USE Surface, ONLY: &
       vegin,        &
       sibwet
  
  USE IOLowLevel, ONLY: &
      ReadVar      , &
      ReadGetNFTGZ
       
  USE InputOutput, ONLY: &
      getsbc, &
      getsbc1D,&
      nfprt , &
      nferr , &
      ifprt 

  USE Options, ONLY: &      
      isimp ,&
      nfcnv0,&
      nfsibt,&
      intgr ,&
      initlz,&		      
      ifalb ,&		      
      ifsst ,&		      
      ifslm ,&
      ifsnw ,&
      sstlag,&
      intsst,&
      fint  ,&
      yrl   ,&      
      monl  ,&
      nftgz0,&
      mxiter,&
      nfsibi
      
  USE GridHistory, ONLY:      &
      IsGridHistoryOn
  IMPLICIT NONE

  ! Gaussian fields: 28 3D , 12 2D and 12 1D

  PRIVATE

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   convc_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   convt_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   convb_in   
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rvisbc_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rvisdc_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rnirbc_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rnirdc_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   avisd_in 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   dswtop_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   gl0_in  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   zorl_in  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   sheleg_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tseam_in   
  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: imask_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rvisb_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rvisd_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rnirb_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   rnird_in  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   gtsea_in 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tc0_in	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg0_in		
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   td0_in		
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: w0_in		
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: capac0_in	 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tcm_in	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tgm_in	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tdm_in	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wm_in	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: capacm_in    
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ppli_in  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ppci_in   
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   var_in
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp1_in 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp2_in  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp3_in  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcpt_in 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   toplv_in 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   botlv_in 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg1_in 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg2_in   
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg3_in    
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   soilm_in
  




  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: sigki 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xvisb 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xvisd 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xnirb 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xnird 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xswtop
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xvisbc
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xvisdc
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xnirbc
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: xnirdc
  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yvisb  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yvisd  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ynirb  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ynird  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yswtop 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yvisbc 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: yvisdc 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ynirbc 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ynirdc 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: cldsav 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ustr
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: vstr
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ssib
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convc
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convt
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convb
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convts
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convcs
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: convbs
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: htrc  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rvisbc	 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rvisdc	 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rnirbc	 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rnirdc	 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: dlwclr
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: uswtpc
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rsclr 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ultclr
  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: avisb 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: avisd
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: anirb 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: anird 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: dswtop
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rs    
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ulwtop  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: gl0   
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: zorl  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: sheleg
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: tseam
  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: htr	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: clr	
  INTEGER(KIND=i8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: imask   
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rvisb 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rvisd 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rnirb 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: rnird 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: dlwbot
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: gtsea
  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tc0	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg0	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   td0	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: w0	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: capac0 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tcm	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tgm	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tdm	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wm	
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: capacm 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ppli   
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ppci
  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   var  

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prct  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcc  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp1 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp2 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcp3 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   prcpt 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   toplv 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   botlv 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   geshem

  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg1  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg2  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tg3  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   soilm 
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   sens
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   evap
  
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: gndvi
  REAL(KIND=r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: ndvim

  INTEGER, PRIVATE :: iMax 
  INTEGER, PRIVATE :: jMax  
  INTEGER, PRIVATE :: klevs
  
  PUBLIC ::  InitFieldsPhyscs 
  PUBLIC ::  InitVariancia
  PUBLIC ::  InitVariancia1D
  PUBLIC ::  InitBoundCond
  PUBLIC ::  InitBoundCond1D
  PUBLIC ::  InitCheckfile
  PUBLIC ::  InitGetsbc
  PUBLIC ::  InitGetsbc1D
CONTAINS

  SUBROUTINE InitFieldsPhyscs(ibMax, kMax, jbMax,iMax_in,jMax_in)
    INTEGER, INTENT(IN) :: ibMax
    INTEGER, INTENT(IN) :: kMax
    INTEGER, INTENT(IN) :: jbMax    
    INTEGER, INTENT(IN) :: iMax_in
    INTEGER, INTENT(IN) :: jMax_in
    iMax=iMax_in
    jMax=jMax_in
    klevs=kMax
    !
    ! map (i,j)
    !
    ALLOCATE(convc_in (iMax,jMax)) 
    convc_in=0.0_r8
    ALLOCATE(convt_in (iMax,jMax)) 
    convt_in=0.0_r8  
    ALLOCATE(convb_in (iMax,jMax))
    convb_in=0.0_r8
    ALLOCATE(rvisbc_in(iMax,jMax))
    rvisbc_in=0.0_r8
    ALLOCATE(rvisdc_in(iMax,jMax))
    rvisdc_in=0.0_r8
    ALLOCATE(rnirbc_in(iMax,jMax))
    rnirbc_in=0.0_r8
    ALLOCATE(rnirdc_in(iMax,jMax))
    rnirdc_in=0.0_r8
    ALLOCATE(avisd_in (iMax,jMax)) 
    avisd_in =0.0_r8
    ALLOCATE(dswtop_in(iMax,jMax)) 
    dswtop_in=0.0_r8  
    ALLOCATE(gl0_in   (iMax,jMax)) 
    gl0_in=0.0_r8 
    ALLOCATE(zorl_in  (iMax,jMax)) 
    zorl_in=0.0_r8 
    ALLOCATE(sheleg_in(iMax,jMax)) 
    sheleg_in=0.0_r8
    ALLOCATE(tseam_in (iMax,jMax)) 
    tseam_in=0.0_r8 
    ALLOCATE(gtsea_in (iMax,jMax)) 	
    gtsea_in =0.0_r8
    ALLOCATE(tc0_in   (iMax,jMax))   
    tc0_in   =0.0_r8
    ALLOCATE(tg0_in   (iMax,jMax))   
    tg0_in   =0.0_r8
    ALLOCATE(td0_in   (iMax,jMax))   
    td0_in   =0.0_r8
    ALLOCATE(w0_in    (iMax,3,jMax)) 
    w0_in    =0.0_r8  
    ALLOCATE(capac0_in(iMax,2,jMax)) 
    capac0_in=0.0_r8
    ALLOCATE(tcm_in   (iMax,jMax))   
    tcm_in   =0.0_r8
    ALLOCATE(tgm_in   (iMax,jMax))   
    tgm_in   =0.0_r8
    ALLOCATE(tdm_in   (iMax,jMax))   
    tdm_in   =0.0_r8 
    ALLOCATE(wm_in    (iMax,3,jMax)) 
    wm_in    =0.0_r8
    ALLOCATE(capacm_in(iMax,2,jMax)) 
    capacm_in=0.0_r8
    ALLOCATE(ppli_in  (iMax,jMax))  
    ppli_in =0.0_r8
    ALLOCATE(ppci_in  (iMax,jMax))   
    ppci_in =0.0_r8
    ALLOCATE(var_in   (iMax,jMax))   
    var_in  =0.0_r8
    ALLOCATE(prcp1_in (iMax,jMax))
    prcp1_in=0.0_r8
    ALLOCATE(prcp2_in (iMax,jMax))  
    prcp2_in=0.0_r8
    ALLOCATE(prcp3_in (iMax,jMax)) 
    prcp3_in=0.0_r8
    ALLOCATE(prcpt_in (iMax,jMax)) 
    prcpt_in=0.0_r8
    ALLOCATE(toplv_in (iMax,jMax)) 
    toplv_in=0.0_r8
    ALLOCATE(botlv_in (iMax,jMax)) 
    botlv_in=0.0_r8
    ALLOCATE(tg1_in   (iMax,jMax))
    tg1_in  =0.0_r8
    ALLOCATE(tg2_in   (iMax,jMax))
    tg2_in  =0.0_r8
    ALLOCATE(tg3_in   (iMax,jMax))
    tg3_in  =0.0_r8
    ALLOCATE(soilm_in (iMax,jMax))
    soilm_in=0.0_r8
    ALLOCATE(imask_in (iMax,jMax)) 	
    imask_in =0
    ALLOCATE(rvisb_in (iMax,jMax))
    rvisb_in=0.0_r8
    ALLOCATE(rnirb_in (iMax,jMax))
    rnirb_in=0.0_r8
    ALLOCATE(rvisd_in (iMax,jMax))
    rvisd_in=0.0_r8
    ALLOCATE(rnird_in (iMax,jMax))
    rnird_in=0.0_r8
    !
    !   map(ib,jb)
    !
    ALLOCATE(sigki (ibMax))
    sigki=0.0_r8
    ALLOCATE(xvisb (ibMax))
    xvisb=0.0_r8
    ALLOCATE(xvisd (ibMax))
    xvisd=0.0_r8
    ALLOCATE(xnirb (ibMax))
    xnirb=0.0_r8
    ALLOCATE(xnird (ibMax))
    xnird=0.0_r8
    ALLOCATE(xswtop(ibMax))
    xswtop=0.0_r8
    ALLOCATE(xvisbc(ibMax))
    xvisbc=0.0_r8
    ALLOCATE(xvisdc(ibMax))
    xvisdc=0.0_r8
    ALLOCATE(xnirbc(ibMax))
    xnirbc=0.0_r8
    ALLOCATE(xnirdc(ibMax))
    xnirdc=0.0_r8
    ALLOCATE(yvisb (ibMax,jbMax)) 
    yvisb =0.0_r8
    ALLOCATE(yvisd (ibMax,jbMax)) 
    yvisd =0.0_r8
    ALLOCATE(ynirb (ibMax,jbMax)) 
    ynirb =0.0_r8
    ALLOCATE(ynird (ibMax,jbMax)) 
    ynird =0.0_r8
    ALLOCATE(yswtop(ibMax,jbMax)) 
    yswtop=0.0_r8
    ALLOCATE(yvisbc(ibMax,jbMax)) 
    yvisbc=0.0_r8
    ALLOCATE(yvisdc(ibMax,jbMax)) 
    yvisdc=0.0_r8
    ALLOCATE(ynirbc(ibMax,jbMax)) 
    ynirbc=0.0_r8
    ALLOCATE(ynirdc(ibMax,jbMax)) 
    ynirdc=0.0_r8
    ALLOCATE(cldsav(ibMax,jbMax)) 
    cldsav=0.0_r8
    ALLOCATE(ustr  (ibMax,jbMax)) 
    ustr=0.0_r8
    ALLOCATE(vstr  (ibMax,jbMax)) 
    vstr=0.0_r8
    ALLOCATE(ssib  (ibMax,jbMax)) 
    ssib=0.0_r8
    ALLOCATE(convc (ibMax,jbMax)) 
    convc=0.0_r8
    ALLOCATE(convt (ibMax,jbMax)) 
    convt=0.0_r8 
    ALLOCATE(convb (ibMax,jbMax)) 
    convb=0.0_r8
    ALLOCATE(convts(ibMax,jbMax))
    convts=0.0_r8
    ALLOCATE(convcs(ibMax,jbMax))
    convcs=0.0_r8
    ALLOCATE(convbs(ibMax,jbMax))
    convbs=0.0_r8
    ALLOCATE(htrc  (ibMax,kMax,jbMax))
    htrc  =0.0_r8 
    ALLOCATE(rvisbc(ibMax,jbMax))
    rvisbc=0.0_r8 
    ALLOCATE(rvisdc(ibMax,jbMax))   	  
    rvisdc=0.0_r8	    	  
    ALLOCATE(rnirbc(ibMax,jbMax))   	  
    rnirbc=0.0_r8	     	      	  
    ALLOCATE(rnirdc(ibMax,jbMax))   	      
    rnirdc=0.0_r8
    ALLOCATE(dlwclr(ibMax,jbMax))
    dlwclr=0.0_r8
    ALLOCATE(uswtpc(ibMax,jbMax))
    uswtpc=0.0_r8
    ALLOCATE(rsclr (ibMax,jbMax))
    rsclr =0.0_r8
    ALLOCATE(ultclr(ibMax,jbMax))
    ultclr=0.0_r8
    ALLOCATE(avisb (ibMax,jbMax)) 
    avisb =0.0_r8     
    ALLOCATE(avisd (ibMax,jbMax)) 
    avisd =0.0_r8
    ALLOCATE(anirb (ibMax,jbMax)) 
    anirb =0.0_r8
    ALLOCATE(anird (ibMax,jbMax)) 
    anird =0.0_r8 
    ALLOCATE(dswtop(ibMax,jbMax)) 
    dswtop=0.0_r8
    ALLOCATE(rs	   (ibMax,jbMax)) 
    rs	  =0.0_r8
    ALLOCATE(ulwtop(ibMax,jbMax)) 
    ulwtop=0.0_r8  
    ALLOCATE(gl0   (ibMax,jbMax)) 
    gl0	  =0.0_r8    
    ALLOCATE(zorl  (ibMax,jbMax)) 
    zorl  =0.0_r8 
    ALLOCATE(sheleg(ibMax,jbMax)) 
    sheleg=0.0_r8
    ALLOCATE(tseam (ibMax,jbMax)) 
    tseam=0.0_r8
    ALLOCATE(htr   (ibMax,kMax,jbMax))  
    htr   =0.0_r8
    ALLOCATE(clr   (ibMax,kMax,jbMax))  
    clr   =0.0_r8
    ALLOCATE(imask (ibMax,jbMax)) 	
    imask =0
    ALLOCATE(rvisb (ibMax,jbMax)) 	
    rvisb =0.0_r8
    ALLOCATE(rvisd (ibMax,jbMax)) 	
    rvisd =0.0_r8
    ALLOCATE(rnirb (ibMax,jbMax)) 	
    rnirb =0.0_r8
    ALLOCATE(rnird (ibMax,jbMax)) 	
    rnird =0.0_r8
    ALLOCATE(dlwbot(ibMax,jbMax)) 	
    dlwbot=0.0_r8
    ALLOCATE(gtsea (ibMax,jbMax)) 	
    gtsea =0.0_r8
    ALLOCATE(tc0   (ibMax,jbMax))   
    tc0   =0.0_r8
    ALLOCATE(tg0   (ibMax,jbMax))   
    tg0   =0.0_r8  
    ALLOCATE(td0   (ibMax,jbMax))   
    td0   =0.0_r8 
    ALLOCATE(w0    (ibMax,3,jbMax)) 
    w0    =0.0_r8 
    ALLOCATE(capac0(ibMax,2,jbMax)) 
    capac0=0.0_r8 
    ALLOCATE(tcm   (ibMax,jbMax))   
    tcm   =0.0_r8  
    ALLOCATE(tgm   (ibMax,jbMax))   
    tgm   =0.0_r8
    ALLOCATE(tdm   (ibMax,jbMax))   
    tdm   =0.0_r8     
    ALLOCATE(wm    (ibMax,3,jbMax)) 
    wm    =0.0_r8
    ALLOCATE(capacm(ibMax,2,jbMax)) 
    capacm=0.0_r8
    ALLOCATE(ppli  (ibMax,jbMax))   
    ppli  =0.0_r8    
    ALLOCATE(ppci  (ibMax,jbMax))   
    ppci  =0.0_r8
    ALLOCATE(var   (ibMax,jbMax))   
    var   =0.0_r8
    ALLOCATE(prct  (ibMax,jbMax))  
    prct   =0.0_r8
    ALLOCATE(prcc  (ibMax,jbMax))  
    prcc    =0.0_r8
    ALLOCATE(prcp1 (ibMax,jbMax))  
    prcp1  =0.0_r8
    ALLOCATE(prcp2 (ibMax,jbMax))  
    prcp2  =0.0_r8
    ALLOCATE(prcp3 (ibMax,jbMax))  
    prcp3  =0.0_r8
    ALLOCATE(prcpt (ibMax,jbMax))  
    prcpt   =0.0_r8
    ALLOCATE(toplv (ibMax,jbMax))  
    toplv   =0.0_r8
    ALLOCATE(botlv (ibMax,jbMax))  
    botlv  =0.0_r8
    ALLOCATE(geshem(ibMax,jbMax))  
    geshem =0.0_r8 
    ALLOCATE(tg1   (ibMax,jbMax))
    tg1	   =0.0_r8
    ALLOCATE(tg2   (ibMax,jbMax))
    tg2    =0.0_r8
    ALLOCATE(tg3   (ibMax,jbMax))
    tg3	   =0.0_r8
    ALLOCATE(soilm (ibMax,jbMax))
    soilm  =0.0_r8
    ALLOCATE(sens  (ibMax,jbMax))
    sens  =0.0_r8
    ALLOCATE(evap  (ibMax,jbMax))
    evap  =0.0_r8


    ALLOCATE(gndvi(ibMax,jbMax))
    gndvi=0.0_r8
    ALLOCATE(ndvim(ibMax,jbMax))
    ndvim=0.0_r8    

 END SUBROUTINE InitFieldsPhyscs

 SUBROUTINE InitVariancia(igwd,nfvar,fNameOrgvar)
    CHARACTER(LEN=*) , INTENT(in   ) :: igwd
    INTEGER          , INTENT(in   ) :: nfvar
    CHARACTER(LEN=*) , INTENT(in   ) :: fNameOrgvar
    OPEN(nfvar, file=TRIM(fNameOrgvar),ACTION="read",FORM="UNFORMATTED")
    IF(igwd.eq.'YES ') THEN
      CALL ReadVar(nfvar,var_in)
      var=var_in
    END IF
 END SUBROUTINE InitVariancia 

 SUBROUTINE InitVariancia1D(igwd,nfvar,fNameOrgvar)
    CHARACTER(LEN=*) , INTENT(in   ) :: igwd
    INTEGER          , INTENT(in   ) :: nfvar
    CHARACTER(LEN=*) , INTENT(in   ) :: fNameOrgvar
    REAL(KIND=r8)                             :: orog
    NAMELIST /MODEL_OROG/orog
    READ (111,MODEL_OROG)
    var=orog
    var_in=orog
 END SUBROUTINE InitVariancia1D 
 
 SUBROUTINE InitBoundCond1D(&
       ibMax,jbMax ,iwrkdm,ityp,ifdy  ,todcld,ids   ,idc   ,ifday ,&
       tod  ,totm ,todsib,radwrk,idate ,idatec,jdt  ,si    ,sl    ,&
       fNameSibmsk,fNameTg3zrl  ,ibMaxPerJB)
       
    INTEGER         , INTENT(IN   ) :: ibMax   
    INTEGER         , INTENT(IN   ) :: jbMax   
    INTEGER         , INTENT(IN   ) :: iwrkdm
    INTEGER         , INTENT(IN   ) :: ityp
    INTEGER         , INTENT(INOUT  ) :: ifdy
    REAL(KIND=r8)            , INTENT(INOUT  ) :: todcld
    INTEGER         , INTENT(INOUT  ) :: ids(:)
    INTEGER         , INTENT(INOUT  ) :: idc(:)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)            , INTENT(IN   ) :: tod
    REAL(KIND=r8)            , INTENT(INOUT  ) :: totm 
    REAL(KIND=r8)            , INTENT(INOUT  ) :: todsib
    REAL(KIND=r8)            , INTENT(OUT  ) :: radwrk(iwrkdm,4)
    INTEGER         , INTENT(IN   ) :: idate(:) 
    INTEGER         , INTENT(IN   ) :: idatec(:)
    INTEGER         , INTENT(IN   ) :: jdt   
    
    REAL(KIND=r8)            , INTENT(IN   ) :: si(:)
    REAL(KIND=r8)            , INTENT(IN   ) :: sl(:)
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibmsk
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameTg3zrl
       
    REAL(KIND=r8)                            :: tice  =271.16e0_r8
    REAL(KIND=r8)                            :: t0    
    REAL(KIND=r8)                            :: sinmax
    INTEGER                         :: j     
    INTEGER                         :: i   
    INTEGER                         :: ncount 
    INTEGER                         :: sibmask
    REAL(KIND=r8)                            :: tsoil1
    REAL(KIND=r8)                            :: tsoil2
    REAL(KIND=r8)                            :: tsoil3
    REAL(KIND=r8)                            :: zorlsoil
    REAL(KIND=r8)                            :: wsib  (ibMax,jbMax)
    REAL(KIND=r8)                            :: zero  =0.0e3_r8
    REAL(KIND=r8)                            :: thousd=1.0e3_r8
    REAL(KIND=r8)            , PARAMETER     :: xl0   =10.0_r8
    NAMELIST /MODEL_VEG/sibmask
    NAMELIST /MODEL_soil/tsoil1 ,tsoil2 ,tsoil3 ,zorlsoil  
    IF(isimp.ne.'YES ') THEN
       convc    =0.0_r8
       convt    =0.0_r8
       convb    =0.0_r8
       prcp1    =0.0_r8
       prcp2    =0.0_r8
       prcp3    =0.0_r8
       prcpt    =0.0_r8
       toplv    =0.0_r8
       botlv    =0.0_r8
       sheleg   =0.0_r8
       sheleg_in=0.0_r8
       
       CALL vegin (si(1) ,sl(1))
       READ (111,MODEL_VEG)              
       imask   =  sibmask
       imask_in=  sibmask            
       !
       !     intgr=2  time integration of surface physical variable is done
       !     by leap-frog implicit scheme. this conseves enegy and h2o.
       !     intgr=1  time integration of surface physical variable is done
       !     by backward implicit scheme.
       !
       intgr =2
       !
       !     initialize sib variables
       !
       IF(ifday.eq.0.and.tod.eq.zero.and.initlz.ge.0) THEN
 	 
          avisd_in =  avisd 
          gtsea_in =  gtsea 
          soilm_in =  soilm 
          sheleg_in=  sheleg
	  
	  CALL getsbc1D (1_i8   ,1_i8     ,avisd_in,gtsea_in,soilm_in,sheleg_in,ifday , &
                        tod  ,idate ,idatec,nfprt ,ifprt ,jdt  ,ifalb    , ifsst, &
			ifslm ,ifsnw ,sstlag,intsst,fint ,tice  ,yrl  ,monl)
	  
	  avisd =  avisd_in 
	  gtsea =  gtsea_in 
	  soilm =  soilm_in 
	  sheleg=  sheleg_in

          READ (111,MODEL_soil)
          
	  tg1 = tsoil1
	  tg2 = tsoil2
	  tg3 = tsoil3
	  zorl= zorlsoil
	  tg1_in =tsoil1
	  tg2_in =tsoil2
	  tg3_in =tsoil3
	  zorl_in=zorlsoil
          t0    =271.17_r8
          sinmax=150.0_r8
          !
          !     use rvisd as temporary for abs(soilm)
          !
          DO j=1,jbMax
             DO i=1,ibMax
                rvisd(i,j)=abs(soilm(i,j))
             END DO
          END DO

          IF(mxiter*ityp.gt.iwrkdm)THEN
             WRITE(nfprt,2435)mxiter,ityp,iwrkdm
             WRITE(nferr,2435)mxiter,ityp,iwrkdm
             STOP 2435
          END IF
         	 
          CALL sibwet(ibMax      ,jbMax      ,rvisd      ,sinmax,imask ,wsib   ,ssib  , &
                      radwrk(1,1),radwrk(1,2),radwrk(1,3),mxiter,nfprt ,ifprt  ,&
		      ibMaxPerJB)
	  WRITE(*,*)'imask,wsib,ssib,var,tg1,tg2,tg3,zorl'
          WRITE(*,*)imask,wsib,ssib,var,tg1,tg2,tg3,zorl
	  ppli=0.0_r8
          ppci=0.0_r8  
          capac0=0.0_r8  
          capacm=0.0_r8  
          !
          !     td0 (deep soil temp) is temporarily defined as tg3
          !
          DO j=1,jbMax
             ncount=0
             DO i=1,ibMax
                gl0(i,j)=xl0
                tseam(i,j)=gtsea(i,j)

                IF(imask(i,j).eq.0) THEN
                   IF(-gtsea(i,j).lt.t0) THEN
                      imask(i,j)=-1
                   END IF
                ELSE
                   ncount=ncount+1
                   w0    (i,1,j)=wsib(i,j)
                   w0    (i,2,j)=wsib(i,j)
                   w0    (i,3,j)=wsib(i,j)
                   td0   (i,  j)=tg3 (i,j)
                   wm    (i,1,j)=wsib(i,j)
                   wm    (i,2,j)=wsib(i,j)
                   wm    (i,3,j)=wsib(i,j)
                   tdm   (i,  j)=tg3 (i,j)
                   tgm   (i,  j)=tg3 (i,j)
                   tcm   (i,  j)=tg3 (i,j)
                   ssib  (i,j  )=0.0_r8
                   IF(soilm(i,j).lt.0.0_r8)ssib(i,j)=wsib(i,j)

                   IF(sheleg(i,j).gt.zero) THEN
                      capac0(i,2,j)=sheleg(i,j)/thousd
                      capacm(i,2,j)=sheleg(i,j)/thousd
                   END IF

                END IF
             END DO
          END DO
       END IF
    END IF
444 FORMAT(' SIB PROGNOSTIC VARIABLES READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)  
2435 FORMAT(' FOR SIBWET MXITER=',I4,' * ITYP=',I3, &
         ' CANNOT EXCEED IWRKDM=',I6)	  
 END SUBROUTINE InitBoundCond1D 

 SUBROUTINE InitBoundCond(&
       ibMax,jbMax ,iwrkdm,ityp,ifdy  ,todcld,ids   ,idc   ,ifday ,&
       tod  ,totm ,todsib,radwrk,idate ,idatec,jdt  ,si    ,sl    ,&
       fNameSibmsk,fNameTg3zrl  ,ibMaxPerJB)
       
    INTEGER         , INTENT(IN   ) :: ibMax   
    INTEGER         , INTENT(IN   ) :: jbMax   
    INTEGER         , INTENT(IN   ) :: iwrkdm
    INTEGER         , INTENT(IN   ) :: ityp
    INTEGER         , INTENT(OUT  ) :: ifdy
    REAL(KIND=r8)            , INTENT(OUT  ) :: todcld
    INTEGER         , INTENT(OUT  ) :: ids(:)
    INTEGER         , INTENT(OUT  ) :: idc(:)
    INTEGER         , INTENT(IN   ) :: ifday
    REAL(KIND=r8)            , INTENT(IN   ) :: tod
    REAL(KIND=r8)            , INTENT(OUT  ) :: totm 
    REAL(KIND=r8)            , INTENT(OUT  ) :: todsib
    REAL(KIND=r8)            , INTENT(OUT  ) :: radwrk(iwrkdm,4)
    INTEGER         , INTENT(IN   ) :: idate(:) 
    INTEGER         , INTENT(IN   ) :: idatec(:)
    INTEGER         , INTENT(IN   ) :: jdt   
    
    REAL(KIND=r8)            , INTENT(IN   ) :: si(:)
    REAL(KIND=r8)            , INTENT(IN   ) :: sl(:)
    INTEGER         , INTENT(IN   ) :: ibMaxPerJB(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameSibmsk
    CHARACTER(LEN=*), INTENT(IN   ) :: fNameTg3zrl
       
    REAL(KIND=r8)                            :: tice  =271.16e0_r8
    REAL(KIND=r8)                            :: t0    
    REAL(KIND=r8)                            :: sinmax
    INTEGER                         :: j     
    INTEGER                         :: i   
    INTEGER                         :: ncount
    REAL(KIND=r8)                            :: wsib  (ibMax,jbMax)
    REAL(KIND=r8)                            :: zero  =0.0e3_r8
    REAL(KIND=r8)                            :: thousd=1.0e3_r8
    REAL(KIND=r8)            , PARAMETER     :: xl0   =10.0_r8
    
    IF(isimp.ne.'YES ') THEN
       IF(nfcnv0.ne.0) THEN
          READ(nfcnv0) ifdy,todcld,ids,idc          
	  READ(nfcnv0) convc,convt,convb,prcp1,prcp2,prcp3, &
                       prcpt,toplv,botlv
          IF(ifday.gt.0.or.tod.gt.0.0_r8)READ(nfcnv0)rvisd,rvisb,rnird, &
                       rnirb,rvisdc,rvisbc,rnirdc,rnirbc,dswtop,totm
          
	  REWIND nfcnv0
	  
       IF(ifprt(4) .ge. 1)WRITE(nfprt,555)ifdy,todcld,ids,idc

       ELSE
          convc=0.0_r8
          convt=0.0_r8
          convb=0.0_r8
          prcp1=0.0_r8
          prcp2=0.0_r8
          prcp3=0.0_r8
          prcpt=0.0_r8
          toplv=0.0_r8
          botlv=0.0_r8
       END IF

       sheleg=0.0_r8
       sheleg_in=0.0_r8
       
       CALL vegin (si(1) ,sl(1))
       OPEN(nfsibt, file=TRIM(fNameSibmsk),ACTION="read",FORM="UNFORMATTED")   
       OPEN(nftgz0, file=TRIM(fNameTg3zrl),ACTION="read",FORM="UNFORMATTED")     
       
       READ  (nfsibt) imask_in
        
       imask=imask_in
       
       REWIND nfsibt 
       ! 
       !     intgr=2  time integration of surface physical variable is done
       !     by leap-frog implicit scheme. this conseves enegy and h2o.
       !     intgr=1  time integration of surface physical variable is done
       !     by backward implicit scheme.
       !
       intgr =2
       !
       !     initialize sib variables
       !
       IF(ifday.eq.0.and.tod.eq.zero.and.initlz.ge.0) THEN
 	 
          avisd_in =avisd 
          gtsea_in =gtsea 
          soilm_in =soilm 
          sheleg_in=sheleg
	  
	  CALL getsbc (iMax ,jMax  ,avisd_in,gtsea_in,soilm_in,sheleg_in,ifday , &
                       tod  ,idate ,idatec,nfprt ,ifprt ,jdt  ,ifalb , &
                       ifsst,ifslm ,ifsnw ,sstlag,intsst,fint ,tice  , &
                       yrl  ,monl)
	       
	  avisd  = avisd_in 
	  gtsea  = gtsea_in 
	  soilm  = soilm_in 
	  sheleg = sheleg_in
	    	  
	  CALL ReadGetNFTGZ(nftgz0,tg1_in,tg2_in,tg3_in,zorl_in)
          
	  tg1  = tg1_in 
	  tg2  = tg2_in 
	  tg3  = tg3_in 
	  zorl = zorl_in

          t0    =271.17_r8
          sinmax=150.0_r8
          !
          !     use rvisd as temporary for abs(soilm)
          !
          DO j=1,jbMax
             DO i=1,ibMaxPerJB(j)
                rvisd(i,j)=abs(soilm(i,j))
             END DO
          END DO

          IF(mxiter*ityp.gt.iwrkdm)THEN
             WRITE(nfprt,2435)mxiter,ityp,iwrkdm
             WRITE(nferr,2435)mxiter,ityp,iwrkdm
             STOP 2435
          END IF
         	 
          CALL sibwet(ibMax ,jbMax ,rvisd ,sinmax   ,imask ,wsib    ,ssib  , &
                      radwrk(1,1),radwrk(1,2),radwrk(1,3),mxiter,nfprt, ifprt,&
		      ibMaxPerJB)

          ppli=0.0_r8
          ppci=0.0_r8  
          capac0=0.0_r8  
          capacm=0.0_r8  
          !
          !     td0 (deep soil temp) is temporarily defined as tg3
          !
          DO j=1,jbMax
             ncount=0
             DO i=1,ibMaxPerJB(j)
                gl0(i,j)=xl0
                tseam(i,j)=gtsea(i,j)

                IF(imask(i,j).eq.0) THEN
                   IF(-gtsea(i,j).lt.t0) THEN
                      imask(i,j)=-1
                   END IF
                ELSE
                   ncount=ncount+1
                   w0    (i,1,j)=wsib(i,j)
                   w0    (i,2,j)=wsib(i,j)
                   w0    (i,3,j)=wsib(i,j)
                   td0   (i,  j)=tg3 (i,j)
                   wm    (i,1,j)=wsib(i,j)
                   wm    (i,2,j)=wsib(i,j)
                   wm    (i,3,j)=wsib(i,j)
                   tdm   (i,  j)=tg3 (i,j)
                   tgm   (i,  j)=tg3 (i,j)
                   tcm   (i,  j)=tg3 (i,j)
                   ssib  (i,j  )=0.0_r8
                   IF(soilm(i,j).lt.0.0_r8)ssib(i,j)=wsib(i,j)

                   IF(sheleg(i,j).gt.zero) THEN
                      capac0(i,2,j)=sheleg(i,j)/thousd
                      capacm(i,2,j)=sheleg(i,j)/thousd
                   END IF

                END IF
             END DO
          END DO

       ELSE

          READ(nfsibi)ifdy,todsib,ids,idc
          READ(nfsibi) td0   ,tdm
          READ(nfsibi) tg0   ,tgm
          READ(nfsibi) tc0   ,tcm
          READ(nfsibi) w0    ,wm
          READ(nfsibi) capac0,capacm
          READ(nfsibi) ppci  ,ppli 
          READ(nfsibi) gl0   ,zorl  ,gtsea ,tseam
	  
          REWIND nfsibi
          	  
	  IF(initlz.lt.0.and.initlz.gt.-3)THEN
	  
	   avisd_in =avisd 
           gtsea_in =gtsea 
           soilm_in =soilm 
           sheleg_in=sheleg
	   
           CALL getsbc(iMax ,jMax ,avisd_in ,gtsea_in ,soilm_in ,sheleg_in,ifday ,tod  , &
                       idate ,idatec,nfprt ,ifprt ,jdt   ,ifalb ,ifsst ,ifslm, &
                       ifsnw ,sstlag,intsst,fint  ,tice  ,yrl   ,monl)
              

	   avisd  =  avisd_in	
	   gtsea  =  gtsea_in	
	   soilm  =  soilm_in	
	   sheleg =  sheleg_in  
	  END IF
	  	     
	  DO j=1,jbMax
             ncount=0
             DO i=1,ibMaxPerJB(j)
                IF(imask(i,j).gt.0)THEN
                   ncount=ncount+1
                   ssib(i,j)=0.0_r8
                   IF(w0(i,1,j).lt.0.0_r8)THEN
                      ssib(i,j)=abs(w0(i,1,j))
                      w0(i,1,j)=abs(w0(i,1,j))
                      w0(i,2,j)=abs(w0(i,2,j))
                      w0(i,3,j)=abs(w0(i,3,j)) 
                      wm(i,1,j)=abs(wm(i,1,j))
                      wm(i,2,j)=abs(wm(i,2,j)) 
                      wm(i,3,j)=abs(wm(i,3,j)) 
                   END IF
                END IF
             END DO
          END DO

          IF(ifprt(5).ge.1)WRITE(nfprt,444) ifdy,todsib,ids,idc

       END IF
    END IF
    
444 FORMAT(' SIB PROGNOSTIC VARIABLES READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)  
2435 FORMAT(' FOR SIBWET MXITER=',I4,' * ITYP=',I3, &
         ' CANNOT EXCEED IWRKDM=',I6)	  
 END SUBROUTINE InitBoundCond 

 SUBROUTINE InitCheckfile(&
       ibMax,jbMax  ,ifdy  ,todcld,ids   ,idc   ,ifday , &
       tod   ,idate ,idatec,jdt   ,todsib,ibMaxPerJB )

    INTEGER, INTENT(IN   ) :: ibMax
    INTEGER, INTENT(IN   ) :: jbMax 
    INTEGER, INTENT(INOUT  ) :: ifdy  
    REAL(KIND=r8)   , INTENT(INOUT  ) :: todcld
    INTEGER, INTENT(INOUT  ) :: ids   (:)    
    INTEGER, INTENT(INOUT  ) :: idc   (:) 
    INTEGER, INTENT(IN   ) :: ifday 
    REAL(KIND=r8)   , INTENT(IN   ) :: tod   
    INTEGER, INTENT(IN   ) :: idate (:) 
    INTEGER, INTENT(IN   ) :: idatec(:)
    INTEGER, INTENT(IN   ) :: jdt   
    REAL(KIND=r8)   , INTENT(INOUT  ) :: todsib
    INTEGER, INTENT(IN   ) :: ibMaxPerJB(:)

    INTEGER                :: j	 
    INTEGER                :: ncount
    INTEGER                :: i	 
    REAL(KIND=r8)                   :: tice  =271.16e0_r8
    !
    !     read cloud dataset for cold start
    !     
       convc=0.0_r8
       convt=0.0_r8
       convb=0.0_r8
       prcp1=0.0_r8
       prcp2=0.0_r8
       prcp3=0.0_r8
       prcpt=0.0_r8
       toplv=0.0_r8
       botlv=0.0_r8

444 FORMAT(' SIB PROGNOSTIC VARIABLES READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)
555 FORMAT(' CLOUD PROGNOSTIC DATA READ IN. AT FORECAST DAY', &
         I8,' TOD ',F8.1/' STARTING',3I3,I5,' CURRENT',3I3,I5)	 

 END SUBROUTINE InitCheckfile

  SUBROUTINE InitGetsbc(ibMax,jbMax,ifday ,tod   ,idate ,idatec,jdt   )
    !
    ! $Author: cptec $
    ! $Date: 2000/01/10 22:27:56 $
    ! $Revision: 1.1 $
    !
    ! getsbc :read surface boundary conditions.
    !
    INTEGER, INTENT(in   ) :: ibMax
    INTEGER, INTENT(in   ) :: jbMax
    INTEGER, INTENT(in   ) :: ifday
    REAL(KIND=r8)   , INTENT(in   ) :: tod
    INTEGER, INTENT(in   ) :: idate (:) 
    INTEGER, INTENT(in   ) :: idatec(:)
    INTEGER, INTENT(in   ) :: jdt

    REAL(KIND=r8)                   :: tice  =271.16e0_r8
    !
    avisd_in = avisd 
    gtsea_in = gtsea 
    soilm_in = soilm 
    sheleg_in= sheleg
    
    CALL getsbc (iMax ,jMax  ,avisd_in,gtsea_in,soilm_in,sheleg_in,ifday , &
    		 tod  ,idate ,idatec,nfprt ,ifprt ,jdt  ,ifalb , &
    		 ifsst,ifslm ,ifsnw ,sstlag,intsst,fint ,tice  , &
    		 yrl  ,monl)
    
    avisd =avisd_in 
    gtsea =gtsea_in 
    soilm =soilm_in 
    sheleg=sheleg_in
 
  END SUBROUTINE InitGetsbc

  SUBROUTINE InitGetsbc1D(ibMax,jbMax,ifday ,tod   ,idate ,idatec,jdt   )
    !
    ! $Author: cptec $
    ! $Date: 2000/01/10 22:27:56 $
    ! $Revision: 1.1 $
    !
    ! getsbc :read surface boundary conditions.
    !
    INTEGER, INTENT(in   ) :: ibMax
    INTEGER, INTENT(in   ) :: jbMax
    INTEGER, INTENT(in   ) :: ifday
    REAL(KIND=r8)   , INTENT(in   ) :: tod
    INTEGER, INTENT(in   ) :: idate (:) 
    INTEGER, INTENT(in   ) :: idatec(:)
    INTEGER, INTENT(in   ) :: jdt

    REAL(KIND=r8)                   :: tice  =271.16e0_r8
    !
    avisd_in =  avisd 
    gtsea_in =  gtsea 
    soilm_in =  soilm 
    sheleg_in=  sheleg
    
    CALL getsbc1D (1_i8    ,1_i8     ,avisd_in,gtsea_in,soilm_in,sheleg_in,ifday , &
    		   tod  ,idate ,idatec  ,nfprt   ,ifprt   ,jdt      ,ifalb , &
    		   ifsst,ifslm ,ifsnw   ,sstlag  ,intsst  ,fint     ,tice  , &
    		   yrl  ,monl)
    
    avisd   =  avisd_in 
    gtsea   =  gtsea_in 
    soilm   =  soilm_in 
    sheleg  =  sheleg_in
 
  END SUBROUTINE InitGetsbc1D
END MODULE FieldsPhyscs
