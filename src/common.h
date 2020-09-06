c
c   $Id: common.h 12 2008-12-11 11:12:34Z ishikawa $
c
      integer writ_db,nwrit_db,writ_3d,nkeisu,
     &  topodt,bnddt,inidt,ndgdt,result,continu,restart,
     &  ntrop,haney,monthly,nvhr,
     &  writ_mn,heatmean,nx10minu,nwr1d,write_1d,result1d,
     &  nwrisop,spinupf,tracerf,tracont,flxmean,vortegf
#ifdef ROKKA
     & ,futidemix
#endif
#if !defined(EXP2007) || !defined(ROKKA)
     & ,resmean
#endif
#ifdef BNDTIDE
     & ,stat_3d
     & ,recfile_2d,step_2d,stat_2d,writ_2d
     & ,recfile_en,step_en,stat_en,writ_en
     & ,recfile_cn,step_cn,stat_cn,writ_cn
     & ,recfile_eq
#endif

      common /jobp/ writ_db,nwrit_db,writ_3d,nkeisu,
     &  topodt,bnddt,inidt,ndgdt,result,continu,restart,
     &  ntrop,haney,monthly,nvhr,
     &  writ_mn,heatmean,nx10minu,nwr1d,write_1d,result1d,
     &  nwrisop,spinupf,tracerf,tracont,flxmean,vortegf
#ifdef ROKKA
     & ,futidemix
#endif
#if !defined(EXP2007) || !defined(ROKKA)
     &     ,resmean
#endif
#ifdef BNDTIDE
     & ,stat_3d
     & ,recfile_2d,step_2d,stat_2d,writ_2d
     & ,recfile_en,step_en,stat_en,writ_en
     & ,recfile_cn,step_cn,stat_cn,writ_cn
     & ,recfile_eq	 
#endif

      integer last,nnmats,matsno,matsn2,nergy,write_3d,nkeis,nkai,
     &  modmtn,modmtnb
     &  ,nsec,cal_year
#if defined(ROKKA) && defined(EXP2007)
     &  ,nsec_scf,nsec_ws,npresm
#endif
#ifdef BNDTIDE
     & ,write_2d,write_en,write_cn
#endif
      common /manage/ last,nnmats,matsno,matsn2,nergy,write_3d,nkeis,
     &  nkai,modmtn,modmtnb
     &  ,nsec,cal_year
#if defined(ROKKA) && defined(EXP2007)
     &  ,nsec_scf(nsf*3),nsec_ws(nws*3),npresm
#endif
#ifdef BNDTIDE
     & ,write_2d,write_en,write_cn
#endif
c
      double precision u,v,t,s,ub,vb,tb,sb,rho,rhoo
      common /totl/ u(im,jm,0:km+1),v(im,jm,0:km+1),
     &  t(im,jm,0:km+1),s(im,jm,0:km+1),
     &  ub(im,jm,0:km+1),vb(im,jm,0:km+1),
     &  tb(im,jm,0:km+1),sb(im,jm,0:km+1),
     &  rho(im,jm,0:km+1),rhoo(im,jm,0:km+1)

      double precision zu,zv,sfund,sfvnd
      common /com1/zu(im,jm),zv(im,jm),sfund(im,jm),sfvnd(im,jm)

      double precision ht,um,vm,hta,uma,vma,htb,umb,
#ifdef BNDTIDE
     &  dump1,dump2,dump3,
#endif
     &  vmb,hu,sfun,sfvn,dtcor,dtcor2,
     &  hrumin,hrvmin,dhru,dhrv,
     &  sfunb,sfvnb,umd,vmd,aindm1,aindm2

      common /trop2/ ht(im,jm),um(im,jm),vm(im,jm),hta(im,jm),
#ifdef BNDTIDE
     &  dump1(im,jm),dump2(im,jm),dump3(im,jm),
#endif
     &  uma(im,jm),vma(im,jm),htb(im,jm),umb(im,jm),
     &  vmb(im,jm),hu(im,jm),sfun(im,jm),sfvn(im,jm),
     &  hrumin(im,jm),hrvmin(im,jm),dhru(im,jm),dhrv(im,jm),
     &  aindm1(im,jm),aindm2(im,jm),sfunb(im,jm),sfvnb(im,jm),
     &  umd(im,jm),vmd(im,jm),dtcor(jm),dtcor2(jm)

      double precision ustar,vstar,hrr,wl,hcl,hcla,hclb,
     &  dzub,dztb,umstar,vmstar,hclu,hclua,hrinv
      common /trop1/ ustar(im,jm,0:km+1),vstar(im,jm,0:km+1),
     &  umstar(im,jm),vmstar(im,jm),hrr(im,jm),wl(im,jm,0:km+1),
     &  dzub(im,jm),dztb(im,jm),
     &  hcl(im,jm),hclb(im,jm),hcla(im,jm),
     &  hclu(im,jm),hclua(im,jm),hrinv(im,jm)
c
c
c
      double precision ua,va,ta,sa,dpdx,dpdy,sfu,sfv,
     1  areat,voltr,volur,rar,dzua,dzumin,dzvmin,xind,yind,
     2  dzu,dzt,tust,sust,tvst,svst,dux,duy,dvx,dvy,volt,
     3  tudf,tvdf,sudf,svdf
      common /com2/ ua(im,jm,0:km+1),va(im,jm,0:km+1),
     1  ta(im,jm,0:km+1),sa(im,jm,0:km+1),dpdx(im,jm,0:km+1),
     2  dpdy(im,jm,0:km+1),sfu(im,jm),sfv(im,jm),
     3  areat(im,jm,0:km+1),voltr(im,jm,0:km+1),
     4  volur(im,jm,0:km+1),rar(im,jm,0:km+1),dzu(im,jm,0:km+1),
     6  dzt(im,jm,0:km+1),tvst(im,jm,0:km+1),
     7  tust(im,jm,0:km+1),sust(im,jm,0:km+1),
     8  svst(im,jm,0:km+1),volt(im,jm,0:km+1),
     9  dux(im,jm,0:km+1),duy(im,jm,0:km+1),
     $  dvx(im,jm,0:km+1),dvy(im,jm,0:km+1),
     $  tudf(im,jm,0:km+1),tvdf(im,jm,0:km+1),
     $  sudf(im,jm,0:km+1),svdf(im,jm,0:km+1),
     $  dzua(im,jm,0:km+1),dzumin(im,jm,0:km+1),dzvmin(im,jm,0:km+1),
     $  xind(im,jm,0:km+1),yind(im,jm,0:km+1)
c
c
c
      double precision td,sd,ud,vd,tked
      common /trid/ td(im,jm,0:km+1),sd(im,jm,0:km+1),
     &   ud(im,jm,0:km+1),vd(im,jm,0:km+1),tked(im,jm,0:km+1)
c
c
      integer ho4,exn,texn
      double precision ex,tex
      common /com4/ ho4(im,jm),exn(im,jm),texn(im,jm)
      common /com4b/ ex(im,jm,0:km+1),tex(im,jm,0:km+1)
c
c
      double precision tm4,sm4,dz,dzr,dzz,dzzr,dp,pd,pm,dep,tvol,
     1  ttvol,cs,csr,tng,sine,cst,cstr,sgn,
     2  cor,ashf,anhf,areaur,areauu
      common /com6/ tm4(0:km+1),sm4(0:km+1),dz(0:km+1),
     &  dzr(0:km+1),dzz(0:km+1),dzzr(0:km+1),
     1  dp(0:km+1),pd(0:km+1),pm(0:km+1),dep(0:km+1),tvol(0:km+1),ttvol,
     2  cs(jm),csr(jm),tng(jm),sine(jm),cst(jm),cstr(jm),
     3  cor(jm),ashf(jm),anhf(jm),areaur(jm),areauu(jm),
     4  sgn(jm)
c
c
c
      double precision pi,omega,radian,radius,grav,dx,dy,dxr,dyr,
     1  dx2r,dy2r,dx4r,dy4r,dxsq,dysq,dxsqr,dysqr,dxdeg,dydeg,slat,
     2  slatu,dtuv,dtts,dttr,dtuv2,dtts2,dttr2,c2dtuv,c2dtts,c2dttr,
     3  hdts,hduv,stan,bsn,bcs,depadp,
     4  ddy,tanfi4,radiur,dxdyr,dydxr,dxddy,gamma,
     5  hupp,vupp,hcnt,vcnt,hcnt2,vcnt2,ddmnar,br,adtts,ahour
     &  ,dttke,c2dttke
      common /com7/ pi,omega,radian,radius,grav,dx,dy,dxr,dyr,dx2r,
     1  dy2r,dx4r,dy4r,dxsq,dysq,dxsqr,dysqr,dxdeg,dydeg,slat,slatu,
     2  dtuv,dtts,dttr,dtuv2,dtts2,dttr2,c2dtuv,c2dtts,c2dttr,hdts,
     3  hduv,stan,bsn,bcs,ddy,tanfi4,radiur,
     4  dxdyr,dydxr,dxddy,gamma(0:km+1),hupp,vupp,hcnt,
     5  vcnt,hcnt2,vcnt2,ddmnar,br,adtts,depadp,ahour
     &  ,dttke,c2dttke

      double precision cxn,cxs,cye,cyw,cne,cse
      common /com8/
     3  cxn(im,jm,0:km+1),cxs(im,jm,0:km+1),cye(im,jm,0:km+1),
     $  cyw(im,jm,0:km+1),cne(im,jm,0:km+1),cse(im,jm,0:km+1)
c
      double precision engexa,engina,engpta,enstro,ttmna,ssmna,ddmna,
     1  tmn,smn,dmn,d2mn,eni,est,stab,spmx,
     2  tterm,sterm,uhad,uvad,uhdf,vhad,vvad,vhdf
      common /enganl/ engexa,engina,engpta,enstro,ttmna,ssmna,ddmna,
     1  tmn(km),smn(km),dmn(0:km+1),d2mn(km),
     &  eni(km),est(km),stab(0:km+1),
     2  spmx(0:km+1),tterm(11),sterm(11),uhad,uvad,uhdf,vhad,vvad,vhdf
c
c
c
c      double precision tini,sini,tref,sref,tref12,sref12,
c     &  wsx,wsy,wflux,fricv3
c         common /bnddt1/ tini(im,jm,0:km+1),sini(im,jm,0:km+1),
c     $  tref(im,jm,0:km+1),sref(im,jm,0:km+1),
c     $  tref12(im,jm,0:km+1,0:imn+1),sref12(im,jm,0:km+1,0:imn+1),
c     &  wsx(im,jm),wsy(im,jm),wflux(im,jm),fricv3(im,jm)
cc
c      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
c     & cloud_d,tot_sol_d,wflux_d,
c     & swr_o,netq_o,netq_i,gref
c      common /bnddt2/ wsx_d(im,jm,nday),wsy_d(im,jm,nday),
c     & sc_wind_d(im,jm,nday),stdev_w_d(im,jm,nday),
c     & dpt_t2m_d(im,jm,nday),
c     & temp_2m_d(im,jm,nday),cloud_d(im,jm,nday),
c     & tot_sol_d(im,jm,nday),wflux_d(im,jm,nday),
c     & swr_o(im,jm),netq_o(im,jm),netq_i(im,jm),
c     & gref(im,jm,0:km+1)
cc
c
      double precision tref,sref,wsx,wsy,
     &  wflux,fricv3,swr_o,netq_o,netq_i
      common /bnddt1/ !!tini(im,jm,0:km+1),sini(im,jm,0:km+1),
     &  tref(im,jm,0:km+1),sref(im,jm,0:km+1),
     &  wsx(im,jm),wsy(im,jm),wflux(im,jm),fricv3(im,jm),
     &  swr_o(im,jm),netq_o(im,jm),netq_i(im,jm)
c
#ifdef GL11M
      double precision tref12,sref12,gref
      common /bnddt2/ tref12(im,jm,0:km+1,0:imn+1),
     &  sref12(im,jm,0:km+1,0:imn+1),
     &  gref(im,jm,0:km+1)
c
      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & cloud_d,tot_sol_d,wflux_d
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),cloud_d(im,jm,nsf),
     & tot_sol_d(im,jm,nsf),wflux_d(im,jm,nsf)
#endif
c
#ifdef EQ100M
      real*4 awind_h,xwind_h,ywind_h,atemp_h,rhmdt_h,
     &  apres_h,preci_h,swrad_h,tref_h,sref_h,
     &  uref_h,vref_h,w1ddata,w0ddata,
     &  wrq,wrw
      common /bnddt3/ awind_h(nhour),xwind_h(nhour),
     &  ywind_h(nhour),atemp_h(nhour),rhmdt_h(nhour),
     &  apres_h(nhour),preci_h(nhour),swrad_h(nhour),
     &  tref_h(nhour,km),sref_h(nhour,km),
     &  uref_h(nhour),vref_h(nhour),w1ddata(km),w0ddata,
     &  wrq(5),wrw(3)
      real*8 gref2,uref,vref,xlat,raintemp,gref,tini,sini
      common /bnddt4/ gref2(im,jm,0:km+1),uref,vref,xlat,raintemp,
     &  gref(im,jm,0:km+1),tini(im,jm,0:km+1),sini(im,jm,0:km+1)
#endif
c
#ifdef PC68M
      real*4 tref12,sref12,grefini,gref,tsref1
      common /bnddt2/ tref12(im,jm,0:km+1,imn),
     &  sref12(im,jm,0:km+1,imn),grefini(jm),
     &  gref(im,jm,0:km+1),tsref1(im,jm,0:km+1)
#ifdef CLIMAT
!======== Climat ===========
#ifdef NCEPSF
!---- NCEP-------
      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & dswrf_d,wflux_d,damflx
#ifdef DLWRF
     & ,dlwrf_d
#else
     & ,cloud_d
#endif
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),
     & dswrf_d(im,jm,nsf),wflux_d(im,jm,nsf),
     & damflx(im,jm)
#ifdef DLWRF
     & ,dlwrf(im,jm,nsf)
#else
     & ,cloud_d(im,jm,nsf)
#endif
!----------------
#else
!-- OMIP
      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & cloud_d,tot_sol_d,wflux_d,damflx
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),cloud_d(im,jm,nsf),
     & tot_sol_d(im,jm,nsf),wflux_d(im,jm,nsf),
     & damflx(im,jm)
! end NCEPSF
#endif
! end CLIMAT
#endif
!end PC68M
#endif
c
c
c
      double precision tkea,tke,tkeb,vdts,vduv,depml,eust,evst,
     &   vdtke,rhodf,rhodh,rhouf,rhouh,rhomix,dleng,vddt,vdds
      common /mlm1/ tkea(im,jm,0:km+1),tke(im,jm,0:km+1),
     &  tkeb(im,jm,0:km+1),vdts(im,jm,0:km+1),vduv(im,jm,0:km+1),
     &  depml(im,jm),eust(im,jm,0:km+1),evst(im,jm,0:km+1),
     &  vdtke(im,jm,0:km+1),rhodf(im,jm,0:km+1),rhodh(im,jm,0:km+1),
     &  rhouf(im,jm,0:km+1),rhouh(im,jm,0:km+1),rhomix(im,jm),
     &  dleng(im,jm,0:km+1),vddt(im,jm,0:km+1),vdds(im,jm,0:km+1)
c
c
c
      double precision ckrmn,rough0,cm0,vdtsmn,vduvmn,calph,
     &   c0,s0,cpr,sigma,tkemin,vdtsmx,drfit,vduvmx,
     &   dleng_e,c0_e,calph_e,pr0,calph_p,calph_b,prandtl
      common /mlm2/ ckrmn,rough0,cm0,vdtsmn,vduvmn,calph,
     &   c0,s0,cpr,sigma,tkemin,vdtsmx,drfit(nmfit),vduvmx,
     &   dleng_e,c0_e,calph_e,pr0,calph_p,calph_b,prandtl
c
c
      double precision aip,adp,ar,agm,drmax,sxymax,hdtsmin,ags,
     &    hdts_btm,hdts_sfc
      common /iso1/ aip,adp,ar,agm,drmax,sxymax,hdtsmin,ags,
     &    hdts_btm,hdts_sfc
c
c
      double precision rho_a,rho_w,rho_i,c_p_a,c_p_w,c_p_i,l_w,l_i,
     &   a_w,b_w,c_w,d_w,a_i,b_i,c_i,d_i,eps_w,cboltz,p_a,
     &   c_la_w,albedo_i,albedo_w,a1_frz,a2_frz,a3_frz,k_i,
     &   hmin,hmax,s_ice,hdice,cd_a_i,cd_a_o,cd_i_o,dtuice,
     &   albedo_melt,c_stab,c_unst,aicemin,dtskin,t_ice_lim,dump_i,
     &   xkai,hminu
      integer modice
      common/bulk/ rho_a,rho_w,rho_i,c_p_a,c_p_w,c_p_i,l_w,l_i,
     &   a_w,b_w,c_w,d_w,a_i,b_i,c_i,d_i,eps_w,cboltz,p_a,
     &   c_la_w,albedo_i,albedo_w,a1_frz,a2_frz,a3_frz,k_i,
     &   hmin,hmax,s_ice,hdice,cd_a_i,cd_a_o,cd_i_o,dtuice,
     &   albedo_melt,c_stab,c_unst,aicemin,dtskin,t_ice_lim,
     &   dump_i,xkai(jm),hminu
      common/bulk_i/ modice
c
c
c
      double precision aicea,aice,aiceb,volicea,volice,voliceb,
     &   uice,vice,hiceu,t_top,ta1toice
      common /ice/ aicea(im,jm),aice(im,jm),aiceb(im,jm),
     &   volicea(im,jm),volice(im,jm),voliceb(im,jm),
     &   uice(im,jm),vice(im,jm),hiceu(im,jm),t_top(im,jm),
     &   ta1toice(im,jm)
c
c
c
      real*8 wmahour
      real*4 wmu,wmv,wmt,wms,wmrho,wmwl,wmhcl,wmsfu,wmsfv,
     &   wmtke,wmvdts,wmmld,wmaice,wmhice,wmuice,wmvice,
     &   wmheat,wmwater,wmvddt,wmvdds,
     &   wmheat_r,wmsal_r
      common /wrm1/ wmahour
      common /wrm2/
     &     wmu(im,jm,0:km+1),wmv(im,jm,0:km+1),wmt(im,jm,0:km+1),
     $     wms(im,jm,0:km+1),wmrho(im,jm,0:km+1),wmwl(im,jm,0:km+1),
     $     wmhcl(im,jm),wmsfu(im,jm),wmsfv(im,jm),
     $     wmtke(im,jm,0:km+1),wmvdts(im,jm,0:km+1),wmmld(im,jm),
     &     wmaice(im,jm),wmhice(im,jm),wmuice(im,jm),wmvice(im,jm),
     &     wmheat(im,jm),wmwater(im,jm),
     &     wmvddt(im,jm,0:km+1),wmvdds(im,jm,0:km+1),
     &     wmheat_r(im,jm),wmsal_r(im,jm)
c
c
c
      real*8 vdtide
      common /tide/ vdtide(im,jm,0:km+1)
c
c
      double precision dratio,expt,exps,vddsmx,dratioc,vddtmn,
     &    vdtsuji
      common /ddifc/ dratio,expt,exps,vddsmx,dratioc,vddtmn,
     &    vdtsuji(0:km+1)
c
c
#ifdef BISMFRIC
      double precision d2ud2x,d2ud2y,d2vd2x,d2vd2y,
     &   bsmags,bsmagw,bsmagn,bsmage,dux2,duy2,dvx2,dvy2,
     &   bsmags_0,bsmagw_0,bsmagn_0,bsmage_0
      common /bihfric/ d2ud2x(im,jm,0:km+1),d2ud2y(im,jm,0:km+1),
     &   d2vd2x(im,jm,0:km+1),d2vd2y(im,jm,0:km+1),
     &   bsmags(im,jm,0:km+1),bsmagw(im,jm,0:km+1),
     &   bsmagn(im,jm,0:km+1),bsmage(im,jm,0:km+1),
     &   dux2(im,jm,0:km+1),duy2(im,jm,0:km+1),
     &   dvx2(im,jm,0:km+1),dvy2(im,jm,0:km+1),
     &   bsmags_0(im,jm),bsmagn_0(im,jm),bsmage_0(im,jm)
     &      ,bsmagw_0(im,jm)
      double precision d2sud2x,d2sud2y,d2svd2x,d2svd2y,delx
      common /bihfric2/ d2sud2x(im,jm),d2sud2y(im,jm),d2svd2x(im,jm),
     &   d2svd2y(im,jm),delx(jm)
c
c      double precision dudxw,dudxe,dudys,dudyn,dvdxw,dvdxe,
c     &   dvdys,dvdyn,defrate_sw,defrate_nw,defrate_se,defrate_ne,
c     &   csmag,hduv_bh
c      common /paramsmag/ dudxw,dudxe,dudys,dudyn,dvdxw,dvdxe,
c     &   dvdys,dvdyn,defrate_sw,defrate_nw,defrate_se,defrate_ne,
c     &   csmag,hduv_bh
c
      double precision csmag,hduv_bh
      common /paramsmag/ csmag,hduv_bh
#endif

!      integer nn,ii,jj,kk,nnset
!      common /debugi/ nn,ii,jj,kk,nnset
!      real*8 xx,xx2
!      common /debugr8/ xx,xx2
c
#ifdef NESTED
#ifdef JP44
      double precision tref12,sref12,gref
      common /bnddt2/ tref12(im,jm,0:km+1,0:imn+1),
     &  sref12(im,jm,0:km+1,0:imn+1),
     &  gref(im,jm,0:km+1)

      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & cloud_d,tot_sol_d,wflux_d
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),cloud_d(im,jm,nsf),
     & tot_sol_d(im,jm,nsf),wflux_d(im,jm,nsf)
#endif

#ifdef JP68M
!++++++++++ JP68M ++++++++++++++++++
      real*4 tref12,sref12,grefini,gref,tsref1
      common /bnddt2/ tref12(im,jm,0:km+1,nrbd),
     &  sref12(im,jm,0:km+1,nrbd),grefini(jm),
     &  gref(im,jm,0:km+1),tsref1(im,jm,0:km+1)

#ifdef CLIMAT
!======== Climat ===========
#ifdef NCEPSF
!---- NCEP-------
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),
     & dswrf_d(im,jm,nsf),wflux_d(im,jm,nsf),
     & damflx(im,jm)
#ifdef DLWRF
     & ,dlwrf(im,jm,nsf)
#else
     & ,cloud_d(im,jm,nsf)
#endif
!----------------
#else
!----- OMIP -----
      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & cloud_d,tot_sol_d,wflux_d,damflx
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),cloud_d(im,jm,nsf),
     & tot_sol_d(im,jm,nsf),wflux_d(im,jm,nsf),
     & damflx(im,jm)
!----------------
#endif
#else
#ifdef NCEPSF
!----- NCEP -----
      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & dswrf_d,wflux_d,damflx
#ifdef DLWRF
     & ,dlwrf
#else
     & ,cloud_d
#endif
      common /bnddt3/ wsx_d(im,jm,0:nsf+1),wsy_d(im,jm,0:nsf+1),
     & sc_wind_d(im,jm,0:nsf+1),stdev_w_d(im,jm,0:nsf+1),
     & dpt_t2m_d(im,jm,0:nsf+1),
     & temp_2m_d(im,jm,0:nsf+1),
     & dswrf_d(im,jm,0:nsf+1),wflux_d(im,jm,0:nsf+1),
     & damflx(im,jm)
#ifdef DLWRF
     & ,dlwrf(im,jm,0:nsf+1)
#else
     & ,cloud_d(im,jm,0:nsf+1)
#endif
!----------------
#endif
!==========================
#endif

      integer contini
      common /unnum/contini
!++++++++++++++++++++++++++++++++++++++++
#endif

#ifdef NWNPAC
!++++++++++++++ NWNPAC ++++++++++++++++++

      double precision gref,tref12,sref12,grefini
      common /bnddt2/ tref12(im,jm,0:km+1,nrbd),
     &  sref12(im,jm,0:km+1,nrbd),
     &  gref(im,jm,0:km+1),grefini(jm)

      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & dswrf_d,wflux_d
#ifdef DLWRF
     & ,dlwrf_d
#else
     & ,cloud_d
#endif
#ifdef CLIMAT
!========== climatology ==========
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),wflux_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),
     & dswrf_d(im,jm,nsf)
#ifdef DLWRF
     & ,dlwrf_d(im,jm,nsf)
#else
     & ,cloud_d(im,jm,nsf)
#endif

#else
      common /bnddt3/ wsx_d(im,jm,0:nsf+1),wsy_d(im,jm,0:nsf+1),
     & sc_wind_d(im,jm,0:nsf+1),stdev_w_d(im,jm,0:nsf+1),
     & dpt_t2m_d(im,jm,0:nsf+1),wflux_d(im,jm,0:nsf+1),
     & temp_2m_d(im,jm,0:nsf+1),
     & dswrf_d(im,jm,0:nsf+1)
#ifdef DLWRF
     & ,dlwrf_d(im,jm,0:nsf+1)
#else
     & ,cloud_d(im,jm,0:nsf+1)
#endif
!==================================
#endif
      integer contini
      common /unnum/contini
!++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

#ifdef ROKKA
!+++++++++++++++ ROkkA +++++++++++++++++++++++++++++++++++++
      double precision gref,tref12,sref12
      common /bnddt2/ tref12(im,jm,0:km+1,nrbd),
     &  sref12(im,jm,0:km+1,nrbd),
     &  gref(im,jm,0:km+1)

c
      real*4 wsx_d,wsy_d,sc_wind_d,stdev_w_d,dpt_t2m_d,temp_2m_d,
     & dswrf_d,wflux_d
#ifdef DLWRF
     & ,dlwrf_d
#else
     & ,cloud_d
#endif

#ifdef CLIMAT
!============ Climatology =================
      common /bnddt3/ wsx_d(im,jm,nsf),wsy_d(im,jm,nsf),
     & sc_wind_d(im,jm,nsf),stdev_w_d(im,jm,nsf),
     & dpt_t2m_d(im,jm,nsf),wflux_d(im,jm,nsf),
     & temp_2m_d(im,jm,nsf),
     & dswrf_d(im,jm,nsf)
#ifdef DLWRF
     & ,dlwrf_d(im,jm,nsf)
#else
     & ,cloud_d(im,jm,nsf)
#endif

#else
#if defined(EXP2006) || defined(EXP2007)
      common /bnddt3/ wsx_d(im,jm,nws*3),wsy_d(im,jm,nws*3),
     & sc_wind_d(im,jm,nsf*3),stdev_w_d(im,jm,nsf*3),
     & dpt_t2m_d(im,jm,nsf*3),wflux_d(im,jm,nsf*3),
     & temp_2m_d(im,jm,nsf*3),
     & dswrf_d(im,jm,nsf*3)
#ifdef DLWRF
     & ,dlwrf_d(im,jm,nsf*3)
#else
     & ,cloud_d(im,jm,nsf*3)
#endif
      integer :: resmean
      common /filenum/resmean(numrm)
#else
      common /bnddt3/ wsx_d(im,jm,nsf+1),wsy_d(im,jm,nsf+1),
     & sc_wind_d(im,jm,nsf+1),stdev_w_d(im,jm,nsf+1),
     & dpt_t2m_d(im,jm,nsf+1),wflux_d(im,jm,nsf+1),
     & temp_2m_d(im,jm,nsf+1),
     & dswrf_d(im,jm,nsf+1)
#ifdef DLWRF
     & ,dlwrf_d(im,jm,nsf+1)
#else
     & ,cloud_d(im,jm,nsf+1)
#endif
#endif
!==================================
#endif
      integer contini
      common /unnum/contini
!++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

      double precision texom,exom,pdom,csom,ahourbd
      integer exnom,texnom,nkaibd
      common /orgndat/exnom(im0,jm0),texnom(im0,jm0),nkaibd(nrbd)

      double precision ubd,vbd,tbd,sbd,tkebd,hclbd,hclubd,dzubd,
     $     dztbd,sfunbd,sfvnbd,umbd,vmbd,xp,yp,zp,sp1,
     $     slon,slon0,slat0,slatu0,dxdeg0,dydeg0,dtbd,
     $     chfht,chfhcl,chft,chfice,chfq,cbf,cay,slonu,slonu0,
     $     uicebd,vicebd,volicebd,aicebd,
     $     dxom,dyom,ddyom,tanfi4om,dyddyom,areatom,ashfom,anhfom,
     $     dxddyom,
     $     sineom,tngom,umom,vmom,dbdvol,rdtr,
     $     wtbd,wtom,
     $     clw,cle,cls,cln
      common /orgdat/
     $     exom(im0,jm0,0:km+1),texom(im0,jm0,0:km+1),
     $     pdom(0:km+1),csom(jm0),areatom(im0,jm0,0:km+1),
     $     ashfom(jm0),anhfom(jm0),sineom(jm0),tngom(jm0),
     $     umom(im0,jm0),vmom(im0,jm0),wtom(im0,jm0)
     $     ,ahourbd(nrbd)

      common /nestdt/
     $     ubd(im,jm,0:km+1),vbd(im,jm,0:km+1),
     $     tbd(im,jm,0:km+1),sbd(im,jm,0:km+1),
     $     tkebd(im,jm,0:km+1),
     $     hclbd(im,jm),hclubd(im,jm),dzubd(im,jm,0:km+1),
     $     sfunbd(im,jm),sfvnbd(im,jm),umbd(im,jm),vmbd(im,jm),
     $     dztbd(im,jm,0:km+1),wtbd(im,jm),
     $     uicebd(im,jm),vicebd(im,jm),volicebd(im,jm),aicebd(im,jm)
      common /intp/ xp(ndpt),yp(ndpt),zp(ndpt)
     $     ,sp1(img+2,jmg+2)

      integer nrstart,nrct,nrunit,mb1,mb2,nrng,nsm,omunt,omtopo,
     $     isumt,isumu,ipt,jpt,ipu,jpu,obc_clim
      common/nestipm/nrstart,nrct,nrunit,
     $   mb1,mb2,nrng,nsm,omunt,omtopo,
     $   isumt(0:km+1),isumu(0:km+1),ipt(ndpt,0:km+1),jpt(ndpt,0:km+1),
     $   ipu(ndpt,0:km+1),jpu(ndpt,0:km+1),obc_clim

      common/nestpm/slon,slon0,slat0,slatu0,dxdeg0,dydeg0,chfht,
     $     chfhcl,chft,chfice,chfq,cbf,cay,slonu,slonu0,dtbd,
     $     dxom,dyom,ddyom,tanfi4om,dyddyom,dbdvol,rdtr,
     $     dxddyom,
     $     clw,cle,cls,cln

      double precision uwbc,uebc,usbc,unbc,vwbc,vebc,vsbc,vnbc,
     $   twbc,tebc,tsbc,tnbc,swbc,sebc,ssbc,snbc,
     $   hclwbc,hclebc,hclsbc,hclnbc,
     $   umwbc,umebc,umsbc,umnbc,vmwbc,vmebc,vmsbc,vmnbc,
     $   tkewbc,tkeebc,tkesbc,tkenbc,
     $   volicewbc,voliceebc,volicesbc,volicenbc,
     $   aicewbc,aiceebc,aicesbc,aicenbc,
     $   uicewbc,vicewbc,uiceebc,viceebc,
     $   uicesbc,vicesbc,uicenbc,vicenbc,
     $   wtwbc,wtebc,wtsbc,wtnbc,
     $   twbt,tebt,tsbt,tnbt,swbt,sebt,ssbt,snbt,
     $   tkewbt,tkeebt,tkesbt,tkenbt,
     $   tkeawbt,tkeaebt,tkeasbt,tkeanbt,
     $   hclwbt,hclebt,hclsbt,hclnbt,hclawbt,hclaebt,hclasbt,hclanbt,
     $   hclbwbt,hclbebt,hclbsbt,hclbnbt,
     $   htwbt,htebt,htsbt,htnbt,htawbt,htaebt,htasbt,htanbt,
     $   htbwbt,htbebt,htbsbt,htbnbt,
     $   umwbt,umebt,umsbt,umnbt,vmwbt,vmebt,vmsbt,vmnbt,
     $   wtbwbt,wtbebt,wtbsbt,wtbnbt,wtwbt,wtebt,wtsbt,wtnbt,
     $   wtawbt,wtaebt,wtasbt,wtanbt,wspwbt,wspebt,wspsbt,wspnbt,
     $   volicewbt,voliceebt,volicesbt,volicenbt,
     $   uicewbt,vicewbt,uiceebt,viceebt,
     $   uicesbt,vicesbt,uicenbt,vicenbt

      common /lbdat1/
     $     uwbc(1:2,jm,0:km+1,nrbd),uebc(1:2,jm,0:km+1,nrbd),
     $     usbc(im,1:2,0:km+1,nrbd),unbc(im,1:2,0:km+1,nrbd)
      common /lbdat2/
     $     vwbc(1:2,jm,0:km+1,nrbd),vebc(1:2,jm,0:km+1,nrbd),
     $     vsbc(im,1:2,0:km+1,nrbd),vnbc(im,1:2,0:km+1,nrbd)
      common /lbdat3/
     $     twbc(1:indg+2,jm,0:km+1,nrbd),
     $     tebc(1:indg+2,jm,0:km+1,nrbd),
     $     tsbc(im,1:jndg+2,0:km+1,nrbd),
     $     tnbc(im,1:jndg+2,0:km+1,nrbd)
      common /lbdat4/
     $     swbc(1:indg+2,jm,0:km+1,nrbd),
     $     sebc(1:indg+2,jm,0:km+1,nrbd),
     $     ssbc(im,1:jndg+2,0:km+1,nrbd),
     $     snbc(im,1:jndg+2,0:km+1,nrbd)
      common /lbdat5/
     $     tkewbc(1:indg+2,jm,0:km+1,nrbd),
     $     tkeebc(1:indg+2,jm,0:km+1,nrbd),
     $     tkesbc(im,1:jndg+2,0:km+1,nrbd),
     $     tkenbc(im,1:jndg+2,0:km+1,nrbd)
      common /lbdat6/
     $     hclwbc(1:indg+2,jm,nrbd),hclebc(1:indg+2,jm,nrbd),
     $     hclsbc(im,1:jndg+2,nrbd),hclnbc(im,1:jndg+2,nrbd),
     $     umwbc(1:2,jm,nrbd),umebc(1:2,jm,nrbd),
     $     umsbc(im,1:2,nrbd),umnbc(im,1:2,nrbd),
     $     vmwbc(1:2,jm,nrbd),vmebc(1:2,jm,nrbd),
     $     vmsbc(im,1:2,nrbd),vmnbc(im,1:2,nrbd),
     $     volicewbc(1:indg+2,jm,nrbd),
     $     voliceebc(1:indg+2,jm,nrbd),
     $     volicesbc(im,1:jndg+2,nrbd),
     $     volicenbc(im,1:jndg+2,nrbd),
     $     aicewbc(1:2,jm,nrbd),aiceebc(1:2,jm,nrbd),
     $     aicesbc(im,1:2,nrbd),aicenbc(im,1:2,nrbd),
     $     uicewbc(1:2,jm,nrbd),uiceebc(1:2,jm,nrbd),
     $     uicesbc(im,1:2,nrbd),uicenbc(im,1:2,nrbd),
     $     vicewbc(1:2,jm,nrbd),viceebc(1:2,jm,nrbd),
     $     vicesbc(im,1:2,nrbd),vicenbc(im,1:2,nrbd),
     $     wtwbc(1:indg+2,jm,nrbd),wtebc(1:indg+2,jm,nrbd),
     $     wtsbc(im,1:jndg+2,nrbd),wtnbc(im,1:jndg+2,nrbd)

      common /lbdat7/
     $     twbt(1:indg+2,jm,0:km+1),tebt(1:indg+2,jm,0:km+1),
     $     tsbt(im,1:jndg+2,0:km+1),tnbt(im,1:jndg+2,0:km+1),
     $     swbt(1:indg+2,jm,0:km+1),sebt(1:indg+2,jm,0:km+1),
     $     ssbt(im,1:jndg+2,0:km+1),snbt(im,1:jndg+2,0:km+1),
     $     tkewbt(1:indg+2,jm,0:km+1),tkeebt(1:indg+2,jm,0:km+1),
     $     tkesbt(im,1:jndg+2,0:km+1),tkenbt(im,1:jndg+2,0:km+1),
     $     tkeawbt(1:indg+2,jm,0:km+1),tkeaebt(1:indg+2,jm,0:km+1),
     $     tkeasbt(im,1:jndg+2,0:km+1),tkeanbt(im,1:jndg+2,0:km+1),
     $     hclwbt(1:indg+2,jm),hclebt(1:indg+2,jm),
     $     hclsbt(im,1:jndg+2),hclnbt(im,1:jndg+2),
     $     hclawbt(1:indg+2,jm),hclaebt(1:indg+2,jm),
     $     hclasbt(im,1:jndg+2),hclanbt(im,1:jndg+2),
     $     hclbwbt(1:indg+2,jm),hclbebt(1:indg+2,jm),
     $     hclbsbt(im,1:jndg+2),hclbnbt(im,1:jndg+2),
     $     htbwbt(1:indg+2,jm),htbebt(1:indg+2,jm),
     $     htbsbt(im,1:jndg+2),htbnbt(im,1:jndg+2),
     $     htwbt(1:indg+2,jm),htebt(1:indg+2,jm),
     $     htsbt(im,1:jndg+2),htnbt(im,1:jndg+2),
     $     htawbt(1:indg+2,jm),htaebt(1:indg+2,jm),
     $     htasbt(im,1:jndg+2),htanbt(im,1:jndg+2),
     $     wtwbt(1:indg+2,jm),wtebt(1:indg+2,jm),
     $     wtsbt(im,1:jndg+2),wtnbt(im,1:jndg+2),
     $     wtawbt(1:indg+2,jm),wtaebt(1:indg+2,jm),
     $     wtasbt(im,1:jndg+2),wtanbt(im,1:jndg+2),
     $     wtbwbt(1:indg+2,jm),wtbebt(1:indg+2,jm),
     $     wtbsbt(im,1:jndg+2),wtbnbt(im,1:jndg+2),
     $     wspwbt(1:indg+2,jm),wspebt(1:indg+2,jm),
     $     wspsbt(im,1:jndg+2),wspnbt(im,1:jndg+2),
     $     volicewbt(1:indg+2,jm),voliceebt(1:indg+2,jm),
     $     volicesbt(im,1:jndg+2),volicenbt(im,1:jndg+2),
     $     uicewbt(1:2,jm),vicewbt(1:2,jm),
     $     uiceebt(1:2,jm),viceebt(1:2,jm),
     $     uicesbt(im,1:2),vicesbt(im,1:2),
     $     uicenbt(im,1:2),vicenbt(im,1:2),
     $     umwbt(1:2,jm,2),umebt(1:2,jm,2),
     $     umsbt(im,1:2,2),umnbt(im,1:2,2),
     $     vmwbt(1:2,jm,2),vmebt(1:2,jm,2),
     $     vmsbt(im,1:2,2),vmnbt(im,1:2,2)
#endif

#ifdef ASSIM
      double precision stdu,stdv,stdt,stds,stdh,avevol,
     &   wgtini,eobst,eobsh,ahourini,exbtt,ectdt,ectds
      integer mdasim,nkaiini,nsecini
      common /obs1/ stdh(im,jm),stdu(im,jm,km),
     &   eobst(im,jm),eobsh(im,jm),
     &   stdv(im,jm,km),stdt(im,jm,km),stds(im,jm,km),
     &  exbtt(im,jm,km),ectdt(im,jm,km),ectds(im,jm,km),
     &  mdasim,avevol,wgtini,nsecini,nkaiini,ahourini

!  tf1,sf1 : to construct data misfit
!  uf2,vf1 : restore for adjoint

      double precision tobs,hobs,xbt_t,ctd_t,ctd_s
	  double precision utrue,vtrue,
     $  ttrue,strue,htrue,qttrue,qstrue,tkeobs,tketrue,tauxtrue,
     $  tauytrue,swtrue,aicetr,volicetr,uicetr,vicetr
      integer nxbt,nctd,i_xbt,j_xbt,i_ctd,j_ctd
      common /obs2/ tobs(im,jm,nobst),hobs(im,jm,nobsh),
     &  xbt_t(kxbt,mx_xbt,nobsp),
     &  ctd_t(kctd,mx_ctd,nobsp),ctd_s(kctd,mx_ctd,nobsp),
     &  nxbt(nobsp),nctd(nobsp),
     &  i_xbt(mx_xbt,nobsp),j_xbt(mx_xbt,nobsp),
     &  i_ctd(mx_ctd,nobsp),j_ctd(mx_ctd,nobsp),
     $  utrue(im,jm,0:km+1),vtrue(im,jm,0:km+1),ttrue(im,jm,0:km+1),
     $  strue(im,jm,0:km+1),htrue(im,jm),qttrue(im,jm),qstrue(im,jm),
     $  tkeobs(im,jm,0:km+1),tketrue(im,jm,0:km+1),tauxtrue(im,jm),
     &  tauytrue(im,jm),swtrue(im,jm),aicetr(im,jm),volicetr(im,jm),
     &  uicetr(im,jm),vicetr(im,jm)

      double precision uinit,vinit,
     $  tinit,sinit,hinit,qtinit,qsinit,tkeinit,tauxinit,
     $  tauyinit,swinit
      common /obs3/ uinit(im,jm,0:km+1),vinit(im,jm,0:km+1),
     &  tinit(im,jm,0:km+1),sinit(im,jm,0:km+1),hinit(im,jm),
     &  qtinit(im,jm),qsinit(im,jm),tkeinit(im,jm,0:km+1),
     &  tauxinit(im,jm),tauyinit(im,jm),swinit(im,jm)

#endif

#ifdef BDOUT
      double precision ahour_od,pd_od,t_od,s_od,hcl_od,
     $     um_od,vm_od,aice_od,
     $     volice_od,uice_od,vice_od,tke_od,ddmnar_od! ,u_od,v_od
      common /wrod1/ ahour_od,ddmnar_od,pd_od(0:km+1),
     $     t_od(im,jm,0:km+1),
     $     s_od(im,jm,0:km+1),hcl_od(im,jm),um_od(im,jm),vm_od(im,im),
     $     aice_od(im,jm),volice_od(im,jm),uice_od(im,jm),
     $     vice_od(im,jm),tke_od(im,jm,0:km+1)
c     $     ,u_od(im,jm,0:km+1),v_od(im,jm,0:km+1)
      integer fuoutdt,nwrod
      common /wrod2/ fuoutdt,nwrod
#endif
#ifdef SGOUT
      double precision ahour_sg,u_sg,v_sg,w_sg,vdts_sg
      common /wrsg1/ahour_sg,u_sg(im,jm,0:km+1),v_sg(im,jm,0:km+1),
     $  w_sg(im,jm,0:km+1),vdts_sg(im,jm,0:km+1)
      integer nwrsg,fuoutsg
      common /wrsg2/nwrsg,fuoutsg
#endif

#ifdef NCIO
      integer ncnt_rs,ncnt_rm
#ifdef BNDTIDE
     $ ,ncnt_2d,ncnt_en,ncnt_eq
#endif

      common /nccom1/ncnt_rs(1),ncnt_rm(1)
#ifdef BNDTIDE
     $ ,ncnt_2d(1),ncnt_en(1),ncnt_eq(1)
#endif

      character(len=70) :: sffile0,sffile1,sffile2
      common /cfn/sffile0,sffile1,sffile2
#endif

#ifdef MASKLVIS
      real(8) :: mask_lvis
      integer :: numlvis
      common/lviscntl/mask_lvis(im,jm),numlvis

#endif

#ifdef BNDTIDE
      integer ut_constituent
      common /ut/ ut_constituent(tcm)
      integer :: resconid,resconmax
      character(len=80) :: rescon, datadir
      common /b/
     &  resconid(ncnfile),rescon(ncnfile),resconmax,datadir
	 
#endif

#ifdef IDEBUG
      double precision eq_ht,eq_uv,eq_id
      common /sg/eq_ht(im,jm),eq_uv(im,jm),eq_id(im,jm)
#endif


