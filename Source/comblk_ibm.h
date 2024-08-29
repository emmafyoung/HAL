C***********************************************************************
C
C       Included file for ibm.F program
C       Defines parameters, variables and common blocks
C
C***********************************************************************
C
C       Array sizes:
C
        integer
     $  im             ,jm             ,kb             ,maxp           ,
     $  mx2            ,nlim           ,maxterms       ,nstage         ,
     $  mx2u           ,curvi
C
C	Define model
C
	logical
     $  polcoms        ,occam          ,pom            ,nemo           ,
     $  circumpolar    ,slevel         ,newname        ,TSopt          ,
     $  Iopt           ,beaching       ,oil
C
C***********************************************************************
C
C       Set size and components of model here:
C
        parameter
     $ (im=380         ,jm=400         ,kb=75)
C
        parameter
     $ (mx2=im*jm      ,maxp=250000    ,nlim=1         ,maxterms=4     ,
     $  nstage=10      ,mx2u=im*jm     ,curvi=2)
C
	parameter
     $ (polcoms=.false.,occam=.false.  ,pom=.false.    ,nemo=.true.    ,
     $  circumpolar=.false.            ,slevel=.false. ,newname=.false.,
     $  TSopt=.false.  ,Iopt=.false.   ,beaching=.false.               ,
     $  oil=.false.)
C
C***********************************************************************
C
C       Constants:
C
        real
     $  Pi             ,Radius_E       ,Deg2Rad
C
        parameter
     $  (Pi=3.14159265  , Radius_E=6.371315E6 ,
     $  Deg2Rad=0.017453293)
C
C-----------------------------------------------------------------------
C
C       Scalars:
C
        real
     $  chst           ,fist           ,dch            ,dfi            ,
     $  Ah             ,Az             ,aell           ,bell           ,
     $  dell           ,stageprogressincrement         ,hatchlength    ,
     $  settlelength   ,swim_base      ,mind           ,maxd           ,
     $  maxdi          ,mind_dia       ,maxd_dia

        integer
     $  Dt             ,nrr            ,ncc            ,ndch           ,
     $  ndfi           ,nz             ,nvl            ,
     $  Npart          ,iseed          ,luini          ,
     $  lurest         ,Dt2            ,velnum         ,itot           ,
     $  true_no_stages ,nzu            ,swimflaga      ,swimflagr

        character
     $  veltyp         ,Named          ,Namet          ,Nameh          ,
     $  mortalityYN    ,growthYN       ,migrationYN    ,diapauseYN     ,
     $  swimming
        character*7
     $  species
C
        common/blkcon/
     $  chst           ,fist           ,dch            ,dfi            ,
     $  Ah             ,Az             ,aell           ,bell           ,
     $  dell           ,stageprogressincrement         ,hatchlength    ,
     $  settlelength   ,swim_base      ,mind           ,maxd           ,
     $  maxdi          ,mind_dia       ,maxd_dia       ,Dt             ,
     $  nrr            ,ncc            ,ndch           ,ndfi           ,
     $  nz             ,nvl            ,Npart          ,
     $  iseed          ,luini          ,lurest         ,Dt2            ,
     $  velnum         ,itot           ,true_no_stages ,nzu            ,
     $  swimflagr      ,swimflaga      ,
     $  veltyp         ,Named          ,Namet          ,Nameh          ,
     $  mortalityYN    ,growthYN       ,migrationYN    ,diapauseYN     ,
     $  swimming       ,species
C
C-----------------------------------------------------------------------
C
C       1-D arrays:
C
        real
     $  Hzeta          ,H              ,Temper         ,cspz           ,
     $  Deppart        ,Deppart_st     ,Pstart_T       ,
     $  Pend_T         ,Part_sw        ,stagelen       ,buoyancy       ,
     $  Salin          ,eggdiam        ,Partwt         ,sealevelrise   ,
     $  Tempav_bio     ,stageprogress  ,bwlarvae       ,hlev           ,
     $  Hzeta_a        ,Hzeta_b        ,lat            ,dlat           ,
     $  Icep           ,VMvel          ,swimu          ,swimv          ,
     $  swim_speed     ,bwegg          ,dia_down       ,dia_up         ,
     $  partd_dia      ,nhours         ,Tempav         ,Salinav        ,
     $  eggd
     
        real*8
     $  Xelev          ,Yelev          ,Xu             ,Yu             ,
     $  Xv             ,Yv             ,Xdeg           ,Ydeg

        integer
     $  iz             ,mz             ,part_st        ,iz1            ,
     $  iz2            ,iu1            ,iu2            ,iv1            ,
     $  iv2            ,biocalcfreq    ,iu             ,mu             ,
     $  mv             ,luinput        ,luoutput       ,dia_flag       ,
     $  bathy_flag

        character
     $  buoytyp        ,stage_spec     ,light          ,VMtyp          ,
     $  diapause_cue
	character*80
     $  velfiles

C
        common/blk1d/
     $  Xelev(Maxp)    ,Yelev(Maxp)    ,Xu(Maxp)       ,Yu(Maxp)       ,
     $  Xv(Maxp)       ,Yv(Maxp)       ,Xdeg(Maxp)     ,Ydeg(Maxp)     ,
     $  Hzeta(Mx2)     ,H(Mx2)         ,Temper(Maxp)   ,Partwt(Maxp)   ,
     $  cspz(Mx2)      ,Deppart(Maxp)  ,Deppart_st(Maxp)               ,
     $  Pstart_T(Maxp) ,Pend_T(Maxp)   ,
     $  Part_sw(Maxp)  ,Salin(Maxp)    ,stagelen(nstage)               ,
     $  buoyancy(nstage)               ,eggdiam(nstage),
     $  sealevelrise(Mx2)              ,Tempav_bio(Maxp)               ,
     $  stageprogress(Maxp)            ,bwlarvae(Maxp) ,hlev(Mx2)      ,
     $  Hzeta_a(Mx2)   ,Hzeta_b(Mx2)   ,lat(jm)        ,dlat(jm)       ,
     $  swimu(Mx2)     ,swimv(Mx2)     ,swim_speed(Maxp)               ,
     $  bwegg(Maxp)    ,mz(Mx2)        ,iz(nlim)       ,part_st(Maxp)  ,
     $  iz1(nlim)      ,iz2(nlim)      ,iu1(nlim)      ,iu2(nlim)      ,
     $  iv1(nlim)      ,iv2(nlim)      ,biocalcfreq(nstage)            ,
     $  iu(nlim)       ,mu(Mx2u)       ,velfiles(20)   ,mv(Mx2u)       ,
     $  luinput(20)    ,luoutput(5)    ,Icep(Maxp)     ,VMvel(nstage)  ,
     $  buoytyp(nstage),stage_spec(nstage)             ,light(Maxp)    ,
     $  VMtyp(nstage)  ,diapause_cue(nstage)           ,
     $  dia_down(nstage)               ,dia_up(nstage) ,dia_flag(Maxp) ,
     $  partd_dia(Maxp),bathy_flag(Maxp)               ,nhours(Maxp)   ,
     $  Tempav(Maxp)   ,Salinav(Maxp)  ,eggd(Maxp)
C
C-----------------------------------------------------------------------
C
C       2-D arrays:
C
        real
     $  lons           ,lats           ,dchcu          ,dfic           
C
        common/blk2d/
     $  lons(im,jm)    ,lats(im,jm)    ,dchcu(im,jm)   ,dfic(im,jm)

C
C-----------------------------------------------------------------------
C
C       3-D arrays: none here - all allocatable
C
C-----------------------------------------------------------------------

















