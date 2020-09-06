!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                   c
!     subroutine set                                                c
!
!   $Id: set.F90 13 2008-12-12 10:53:49Z ishikawa $
!                                                                   c
!     This subroutine sets up fundamental model parameters          c
!     including horizontal and vertical resolutions.                c
!                                                                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set


      use param
      use mod_mpi
      use mod_time
#ifdef NCIO
      use netcdf
#endif
#ifdef ASSIM
      use mod_assim
#endif
#ifdef BNDTIDE
      use mod_nao
      use mod_dump
#endif
      implicit none

#include "common.h"

      integer lateq,ir_buf(im,jmg)
      integer i,j,k,n,jlc,idep,idpst,nnset,ii,jj,kk,nfj,m,count1,count2

      real*8 t_frz
      real*8 cvdtide
#ifdef ASSIM
      real*8 dum_ar(2)
      real*8 sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      integer ierr,ipxx,ipyy,ii2,jj2
#else
!  set jobp namelist

!     logocal nfirst : default .true. -> see mod_time.F90
!      logical set_start -> see mod_time.F90
      integer step_3d_year,step_3d_month,step_3d_day,
     &    step_3d_hour,step_3d_min,step_3d_sec,
     & stat_3d_year,stat_3d_month,stat_3d_day,
     &   stat_3d_hour,stat_3d_min,stat_3d_sec,
     & writ_3d_year,writ_3d_month,writ_3d_day,
     &   writ_3d_hour,writ_3d_min,writ_3d_sec,
     & writ_db_year,writ_db_month,writ_db_day,
     &   writ_db_hour,writ_db_min,writ_db_sec,
     & writ_mn_year,writ_mn_month,writ_mn_day,
     &   writ_mn_hour,writ_mn_min,writ_mn_sec,
     & nkeisu_year,nkeisu_month,nkeisu_day,
     &   nkeisu_hour,nkeisu_min,nkeisu_sec,
     & nstart_year,nstart_month,nstart_day,
     &   nstart_hour,nstart_min
#ifdef BNDTIDE
     & ,writ_2d_year,writ_2d_month,writ_2d_day
     &   ,writ_2d_hour,writ_2d_min,writ_2d_sec
     & ,writ_en_year,writ_en_month,writ_en_day
     &   ,writ_en_hour,writ_en_min,writ_en_sec
     & ,writ_cn_year,writ_cn_month,writ_cn_day
     &   ,writ_cn_hour,writ_cn_min,writ_cn_sec
     & ,step_2d_year,step_2d_month,step_2d_day
     &   ,step_2d_hour,step_2d_min,step_2d_sec
     & ,stat_2d_year,stat_2d_month,stat_2d_day
     &   ,stat_2d_hour,stat_2d_min,stat_2d_sec
     & ,step_en_year,step_en_month,step_en_day
     &   ,step_en_hour,step_en_min,step_en_sec
     & ,stat_en_year,stat_en_month,stat_en_day
     &   ,stat_en_hour,stat_en_min,stat_en_sec
     & ,step_cn_year,step_cn_month,step_cn_day
     &   ,step_cn_hour,step_cn_min,step_cn_sec
     & ,stat_cn_year,stat_cn_month,stat_cn_day
     &   ,stat_cn_hour,stat_cn_min,stat_cn_sec
#endif
      real(8) start_sec

      namelist /jobp/ nfirst,
     & step_3d_year,step_3d_month,step_3d_day,
     &   step_3d_hour,step_3d_min,step_3d_sec,
     & stat_3d_year,stat_3d_month,stat_3d_day,
     &   stat_3d_hour,stat_3d_min,stat_3d_sec,
     & writ_3d_year,writ_3d_month,writ_3d_day,
     &   writ_3d_hour,writ_3d_min,writ_3d_sec,
     & writ_db_year,writ_db_month,writ_db_day,
     &   writ_db_hour,writ_db_min,writ_db_sec,
     & writ_mn_year,writ_mn_month,writ_mn_day,
     &   writ_mn_hour,writ_mn_min,writ_mn_sec,
     & nkeisu_year,nkeisu_month,nkeisu_day,
     &   nkeisu_hour,nkeisu_min,nkeisu_sec,
     & set_start,nstart_year,nstart_month,nstart_day,
     &   nstart_hour,nstart_min,start_sec
#ifdef BNDTIDE
     & ,writ_2d_year,writ_2d_month,writ_2d_day
     &   ,writ_2d_hour,writ_2d_min,writ_2d_sec
     & ,step_2d_year,step_2d_month,step_2d_day
     &   ,step_2d_hour,step_2d_min,step_2d_sec
     & ,stat_2d_year,stat_2d_month,stat_2d_day
     &   ,stat_2d_hour,stat_2d_min,stat_2d_sec
     & ,step_en_year,step_en_month,step_en_day
     &   ,step_en_hour,step_en_min,step_en_sec
     & ,stat_en_year,stat_en_month,stat_en_day
     &   ,stat_en_hour,stat_en_min,stat_en_sec
     & ,writ_en_year,writ_en_month,writ_en_day
     &   ,writ_en_hour,writ_en_min,writ_en_sec
     & ,step_cn_year,step_cn_month,step_cn_day
     &   ,step_cn_hour,step_cn_min,step_cn_sec
     & ,stat_cn_year,stat_cn_month,stat_cn_day
     &   ,stat_cn_hour,stat_cn_min,stat_cn_sec
     & ,writ_cn_year,writ_cn_month,writ_cn_day
     &   ,writ_cn_hour,writ_cn_min,writ_cn_sec
#endif

#ifdef BNDTIDE
      integer ut_2n2,ut_2q1,ut_alp2,ut_bet2,ut_dlt2,ut_eps2,ut_et2,
     &   ut_et2s,ut_gam2,ut_j1,ut_j1s,ut_k1,ut_k1s,ut_k1s2,ut_k2,
     &   ut_k2s,ut_kai1,ut_l2,ut_lmd2,ut_m1,ut_m1s,ut_m2,ut_m2s,
     &   ut_mf,ut_mfs,ut_mm,ut_mms1,ut_mms2,ut_msf,ut_msm,ut_mtm,
     &   ut_mu2,ut_mu2s,ut_n2,ut_n2s,ut_nu2,ut_nu2s,ut_o1,ut_o1s,
     &   ut_oo1,ut_oo1s,ut_p1,ut_p1s,ut_phi1,ut_pi1,ut_psi1,ut_q1,
     &   ut_q1s,ut_r2,ut_rho1,ut_s2,ut_s2s,ut_sa,ut_sgm1,ut_ssa,
     &   ut_ssas,ut_t2,ut_tau1,ut_the1,ut_zet2,ut_mtms
      namelist/ut/ ut_2n2,ut_2q1,ut_alp2,ut_bet2,ut_dlt2,ut_eps2,ut_et2,
     &   ut_et2s,ut_gam2,ut_j1,ut_j1s,ut_k1,ut_k1s,ut_k1s2,ut_k2,
     &   ut_k2s,ut_kai1,ut_l2,ut_lmd2,ut_m1,ut_m1s,ut_m2,ut_m2s,
     &   ut_mf,ut_mfs,ut_mm,ut_mms1,ut_mms2,ut_msf,ut_msm,ut_mtm,
     &   ut_mu2,ut_mu2s,ut_n2,ut_n2s,ut_nu2,ut_nu2s,ut_o1,ut_o1s,
     &   ut_oo1,ut_oo1s,ut_p1,ut_p1s,ut_phi1,ut_pi1,ut_psi1,ut_q1,
     &   ut_q1s,ut_r2,ut_rho1,ut_s2,ut_s2s,ut_sa,ut_sgm1,ut_ssa,
     &   ut_ssas,ut_t2,ut_tau1,ut_the1,ut_zet2,ut_mtms
#endif

#endif
#ifdef NCIO
      integer :: ncstat,nsfr,dimid_day
#ifdef ASSIM
      character(len=80) :: rsdir,rscase,incase,workdir
      namelist /filename/  rsdir,rscase,incase,workdir
      character(len=80) :: outfile
#else

#if defined(ROKKA) && defined(EXP2007)
!     character(len=80) :: rsfile,rmfile(numrm),flxmfile,
      character(len=80) :: rsfile,flxmfile,
     &     wsfile0,wsfile1,wsfile2
#else
!     character(len=80) :: rsfile,rmfile,flxmfile
      character(len=80) :: rsfile,flxmfile
#endif
#ifdef BNDTIDE
      character(len=80) :: filename_2d,filename_en,
     &                     filepath_3d,filename_cn(ncnfile),filepath_cn,
     &                     filename_eq
      character(len=2)  :: constr
#endif

!     namelist /filename/ rsfile,rmfile,sffile0,sffile1,sffile2,flxmfile
      namelist /filename/ rsfile,sffile0,sffile1,sffile2,flxmfile
#ifdef BNDTIDE
     &      ,filepath_3d,filename_2d,filename_en,filepath_cn
     &      ,filename_eq
     &      ,datadir
#endif
#if defined(ROKKA) && (defined(EXP2006) || defined(EXP2007))
     &     ,wsfile0,wsfile1,wsfile2
      integer :: l,fdid,nn
#endif
#endif
#endif

#ifndef EVAPO
      real(4) :: evapo_d(im,jm,0:nsf+1)
#endif

!    set job parameters
!
!     dtuv: time of one step in second for momentum equations.
!     dtts: time of one step in second for tracer equations.
!     dttr: time of one step in second for barotropic equations.

#ifdef GL11M
      dtuv = 60.d0*80.d0
      dtts = 60.d0*80.d0
      dttr = 60.d0
      dttke= 60.d0*10.d0
#endif
#ifdef EQ100M
      dtuv=10.d0
      dtts=10.d0
      dttr=2.d0
      dttke=10.d0
#endif
#ifdef PC68M
!      dtuv=60.d0*2.d0
!      dtts=60.d0*2.d0
!      dttr=3.d0
!      dtuv=60.d0*6.d0
!      dtts=60.d0*6.d0
!      dttr=6.d0
       dtuv=60.d0*10.d0
      dtts=60.d0*10.d0
      dttr=10.d0
      dttke=30.d0
#endif
#ifdef JP44
!      dtuv = 60.d0*18.d0
!      dtts = 60.d0*18.d0
!      dttr = 8.d0
      dttke=60.d0*2.d0
      dtuv = 60.d0*18.d0
      dtts = 60.d0*18.d0
      dttr = 8.d0
!      dttke=30.d0
#endif
#ifdef NWNPAC
!----------- initial -----------------
!      dtuv = 60.d0
!      dtts = 60.d0
!      dttr = 1.d0
!      dttke= 60.d0
      dtuv = 90.d0
      dtts = 90.d0
      dttr = 2.d0
      dttke= 90.d0
!      dtuv = 120.d0
!      dtts = 120.d0
!      dttr = 2.d0
!      dttke= 120.d0
c      dtuv = 150.d0
c      dttke= 30.d0
c      dtts = 150.d0
c      dttr = 2.d0
#endif
#ifdef JP68M
!cc   initial:: dtuv/ts: 4min ok 6min ng, dttr: 6ng? 8sec ng.
!      dtuv=60.d0*3.d0
!      dtts=60.d0*3.d0
!      dttke=60.d0*3.d0
!      dttr=4.d0
!cc   after 2 months:: dtuv/ts: 8min ok 10 ng, dttr: 8sec ok 10 ng
      dtuv=60.d0*4.d0
      dtts=60.d0*4.d0
!      dtuv=60.d0*3.d0
!      dtts=60.d0*3.d0
!      dtuv=60.d0*2.d0
!      dtts=60.d0*2.d0
      dttr=4.d0
!      dttr=10.d0
!      dttke=60.d0*8.d0
!      dttr=3.d0
      dttke=60.d0
#endif
#ifdef ROKKA
!      dtuv=1.d0
!      dtts=1.d0
!      dttr=1.d-1
!      dttke=1.d0
!------------ initial ------------
!      dtuv=15.d0
!      dtts=15.d0
!      dttr=5.d-1
!      dttke=15.d0
       dtuv=20.d0
       dtts=20.d0
       dttr=5.d-1
       dttke=20.d0
#ifdef EXPTIDE
!      dtuv=5.d0
!      dtts=5.d0
!      dttr=1.d-1
!      dttke=5.d0 
!      dtuv=20.d0
!      dtts=20.d0
!      dttr=5.d-1
!      dttke=20.d0 
       dtuv=60.d0
       dtts=60.d0
       dttr=5.d-1
       dttke=60.d0
#endif
       dtuv=60.d0
       dtts=60.d0
       dttr=5.d-1
       dttke=60.d0
#endif

      diff_dtuv=timediff(dtuv)
      diff_dttr=timediff(dttr)


!   see readme.txt
!
!     nfirst: .true. (DEFAULT): the present job starts from the initial
!             t/s and state at rest.
!             .false: the present job starts from the last
!             state in the file specified by instrument number 'start'.
! Following param are contains (int) year, month,day,hour,min,sec
!   default values are 0 for all
!     step_3d : number of time steps that the present job should calculate
!     writ_db  : mean values are diagnosed once every 'writ_db' steps.
!             every 'nwrit_db' steps.
!     writ_3d : a set of data is written on the file once every 'writ_3d'
!             steps.
!     writ_2d :a set of tidal data is written on the file once every
!             'writ_2d' steps.
!     nkeisu: coefficients in the equation of state including pressure
!             are recalculated every 'nkeisu' steps.
!    start: must be set if nfirst==.true.
!      (logical) set_start : .true. set start date&time

! set default value
      nfirst=.true.
      step_3d_year=0 ; step_3d_month=0 ; step_3d_day=0
        step_3d_hour=0 ; step_3d_min=0 ; step_3d_sec=0
      stat_3d_year=0 ; stat_3d_month=0 ; stat_3d_day=0
        stat_3d_hour=0 ; stat_3d_min=0 ; stat_3d_sec=0
      writ_3d_year=0 ; writ_3d_month=0 ; writ_3d_day=0
        writ_3d_hour=0 ; writ_3d_min=0 ; writ_3d_sec=0
      writ_db_year=0 ; writ_db_month=0 ; writ_db_day=0
        writ_db_hour=0 ; writ_db_min=0 ; writ_db_sec=0
      writ_mn_year=0 ; writ_mn_month=0 ; writ_mn_day=0
        writ_mn_hour=0 ; writ_mn_min=0 ; writ_mn_sec=0
      nkeisu_year=0 ; nkeisu_month=0 ; nkeisu_day=0
        nkeisu_hour=0 ; nkeisu_min=0 ; nkeisu_sec=0
!***
#ifdef BNDTIDE
      writ_2d_min=0;
      step_2d_year=0 ; step_2d_month=0 ; step_2d_day=0
        step_2d_hour=0 ; step_2d_min=0 ; step_2d_sec=0
      stat_2d_year=0 ; stat_2d_month=0 ; stat_2d_day=0
        stat_2d_hour=0 ; stat_2d_min=0 ; stat_2d_sec=0
      step_en_year=0 ; step_en_month=0 ; step_en_day=0
        step_en_hour=0 ; step_en_min=0 ; step_en_sec=0
      stat_en_year=0 ; stat_en_month=0 ; stat_en_day=0
        stat_en_hour=0 ; stat_en_min=0 ; stat_en_sec=0
      ut_2n2=0
      ut_2q1=0
      ut_alp2=0
      ut_bet2=0
      ut_dlt2=0
      ut_eps2=0
      ut_et2=0
      ut_et2s=0
      ut_gam2=0
      ut_j1=0
      ut_j1s=0
      ut_k1=0
      ut_k1s=0
      ut_k1s2=0
      ut_k2=0
      ut_k2s=0
      ut_kai1=0
      ut_l2=0
      ut_lmd2=0
      ut_m1=0
      ut_m1s=0
      ut_m2=0
      ut_m2s=0
      ut_mf=0
      ut_mfs=0
      ut_mm=0
      ut_mms1=0
      ut_mms2=0
      ut_msf=0
      ut_msm=0
      ut_mtm=0
      ut_mtms=0
      ut_mu2=0
      ut_mu2s=0
      ut_n2=0
      ut_n2s=0
      ut_nu2=0
      ut_nu2s=0
      ut_o1=0
      ut_o1s=0
      ut_oo1=0
      ut_oo1s=0
      ut_p1=0
      ut_p1s=0
      ut_phi1=0
      ut_pi1=0
      ut_psi1=0
      ut_q1=0
      ut_q1s=0
      ut_r2=0
      ut_rho1=0
      ut_s2=0
      ut_s2s=0
      ut_sa=0
      ut_sgm1=0
      ut_ssa=0
      ut_ssas=0
      ut_t2=0
      ut_tau1=0
      ut_the1=0
      ut_zet2=0

#endif

      set_start = .false.
      nstart_year=0 ; nstart_month=0 ; nstart_day=0
      nstart_hour=0 ; nstart_min=0 ; start_sec=0.

      nxday=60*60*24/nint(dtuv)
      nxyear=nday*60*60*24/nint(dtuv)
      nxmonth=nxyear/12
      nxhour=60*60/nint(dtuv)
      nx10minu=60*10/nint(dtuv)
      nxmin = 60/nint(dtuv)
      nxsec = 1./dtuv

#ifndef ASSIM
      if(ip==imaster) then
        read(99,nml=jobp)

! tentative before using mod_datetime
      step_3d = nxyear*step_3d_year+nxmonth*step_3d_month
     &     +nxday*step_3d_day+nxhour*step_3d_hour+nxmin*step_3d_min
     &     +int(nxsec*step_3d_sec)
      writ_3d = nxyear*writ_3d_year+nxmonth*writ_3d_month
     &     +nxday*writ_3d_day+nxhour*writ_3d_hour+nxmin*writ_3d_min
     &     +int(nxsec*writ_3d_sec)
      stat_3d = nxyear*stat_3d_year+nxmonth*stat_3d_month
     &     +nxday*stat_3d_day+nxhour*stat_3d_hour+nxmin*stat_3d_min
     &     +int(nxsec*stat_3d_sec)
      writ_db = nxyear*writ_db_year+nxmonth*writ_db_month
     &     +nxday*writ_db_day+nxhour*writ_db_hour+nxmin*writ_db_min
     &     +int(nxsec*writ_db_sec)
      writ_mn=nxyear*writ_mn_year+nxmonth*writ_mn_month
     &     +nxday*writ_mn_day+nxhour*writ_mn_hour+nxmin*writ_mn_min
     &     +int(nxsec*writ_mn_sec)
      nkeisu= nxyear*nkeisu_year+nxmonth*nkeisu_month+nxday*nkeisu_day
     &     +nxhour*nkeisu_hour

!***
#ifdef BNDTIDE

      read(99,nml=ut)
      writ_2d  = nxyear*writ_2d_year+nxmonth*writ_2d_month
     &     +nxday*writ_2d_day+nxhour*writ_2d_hour+nxmin*writ_2d_min
     &     +int(nxsec*writ_2d_sec)
      step_2d  = nxyear*step_2d_year+nxmonth*step_2d_month
     &     +nxday*step_2d_day+nxhour*step_2d_hour+nxmin*step_2d_min
     &     +int(nxsec*step_2d_sec)
      stat_2d  = nxyear*stat_2d_year+nxmonth*stat_2d_month
     &     +nxday*stat_2d_day+nxhour*stat_2d_hour+nxmin*stat_2d_min
     &     +int(nxsec*stat_2d_sec)
      writ_en  = nxyear*writ_en_year+nxmonth*writ_en_month
     &     +nxday*writ_en_day+nxhour*writ_en_hour+nxmin*writ_en_min
     &     +int(nxsec*writ_en_sec)
      step_en  = nxyear*step_en_year+nxmonth*step_en_month
     &     +nxday*step_en_day+nxhour*step_en_hour+nxmin*step_en_min
     &     +int(nxsec*step_en_sec)
      stat_en  = nxyear*stat_en_year+nxmonth*stat_en_month
     &     +nxday*stat_en_day+nxhour*stat_en_hour+nxmin*stat_en_min
     &     +int(nxsec*stat_en_sec)
      writ_cn  = nxyear*writ_cn_year+nxmonth*writ_cn_month
     &     +nxday*writ_cn_day+nxhour*writ_cn_hour+nxmin*writ_cn_min
     &     +int(nxsec*writ_cn_sec)
      step_cn  = nxyear*step_cn_year+nxmonth*step_cn_month
     &     +nxday*step_cn_day+nxhour*step_cn_hour+nxmin*step_cn_min
     &     +int(nxsec*step_cn_sec)
      stat_cn  = nxyear*stat_cn_year+nxmonth*stat_cn_month
     &     +nxday*stat_cn_day+nxhour*stat_cn_hour+nxmin*stat_cn_min
     &     +int(nxsec*stat_cn_sec)
      read(99,nml=filename)
      resconid(:) = 0
      rescon(:) = ''
      if((writ_cn.gt.0).and.(step_cn.gt.stat_cn)) then
         resconmax = int((step_cn-stat_cn)/writ_cn)
         write(*,*)'resconmax',resconmax
         do n=1,resconmax
            resconid(n)=100+n
            call i2s(n,2,constr)
            rescon(n) = trim(filepath_cn)//constr//'.dat'
            write(*,*)'rescon',rescon(n)
         enddo
      else
         resconmax = 0
      endif

      ut_constituent(1)=ut_q1   ! 16 major short tides
      ut_constituent(2)=ut_o1
      ut_constituent(3)=ut_m1
      ut_constituent(4)=ut_p1
      ut_constituent(5)=ut_k1
      ut_constituent(6)=ut_j1
      ut_constituent(7)=ut_oo1
      ut_constituent(8)=ut_2n2
      ut_constituent(9)=ut_mu2
      ut_constituent(10)=ut_n2
      ut_constituent(11)=ut_nu2
      ut_constituent(12)=ut_m2
      ut_constituent(13)=ut_l2
      ut_constituent(14)=ut_t2
      ut_constituent(15)=ut_s2
      ut_constituent(16)=ut_k2

      ut_constituent(17)=ut_2q1   ! 33 minor short tides
      ut_constituent(18)=ut_sgm1
      ut_constituent(19)=ut_q1s
      ut_constituent(20)=ut_rho1
      ut_constituent(21)=ut_o1s
      ut_constituent(22)=ut_tau1
      ut_constituent(23)=ut_m1s
      ut_constituent(24)=ut_kai1
      ut_constituent(25)=ut_pi1
      ut_constituent(26)=ut_p1s
      ut_constituent(27)=ut_k1s
      ut_constituent(28)=ut_k1s2
      ut_constituent(29)=ut_psi1
      ut_constituent(30)=ut_phi1
      ut_constituent(31)=ut_the1
      ut_constituent(32)=ut_j1s
      ut_constituent(33)=ut_oo1s
      ut_constituent(34)=ut_eps2
      ut_constituent(35)=ut_mu2s
      ut_constituent(36)=ut_n2s
      ut_constituent(37)=ut_nu2s
      ut_constituent(38)=ut_gam2
      ut_constituent(39)=ut_alp2
      ut_constituent(40)=ut_m2s
      ut_constituent(41)=ut_bet2
      ut_constituent(42)=ut_dlt2
      ut_constituent(43)=ut_lmd2
      ut_constituent(44)=ut_s2s
      ut_constituent(45)=ut_r2
      ut_constituent(46)=ut_k2s
      ut_constituent(47)=ut_zet2
      ut_constituent(48)=ut_et2
      ut_constituent(49)=ut_et2s

      ut_constituent(50)=ut_mtm   ! 7 major long tides
      ut_constituent(51)=ut_mf
      ut_constituent(52)=ut_msf
      ut_constituent(53)=ut_mm
      ut_constituent(54)=ut_msm
      ut_constituent(55)=ut_ssa
      ut_constituent(56)=ut_sa 

      ut_constituent(57)=ut_mtms  ! 5 minor long tides
      ut_constituent(58)=ut_mfs
      ut_constituent(59)=ut_mms1
      ut_constituent(60)=ut_mms2
      ut_constituent(61)=ut_ssas


#endif

      endif

	  
      if( set_start ) then
        date_start = init_date1(nstart_year,nstart_month,nstart_day,
     &     nstart_hour,nstart_min,start_sec)
      write(*,*)ip,'start date :',date_start
      endif

	  
      call bcast_logical1(nfirst)
      call bcast_int1(step_3d,1)
      call bcast_int1(stat_3d,1)
      call bcast_int1(writ_3d,1)
      call bcast_int1(writ_db,1)
      call bcast_int1(step_2d,1)
      call bcast_int1(stat_2d,1)
      call bcast_int1(writ_2d,1)
      call bcast_int1(step_en,1)
      call bcast_int1(stat_en,1)
      call bcast_int1(writ_en,1)
      call bcast_int1(step_cn,1)
      call bcast_int1(stat_cn,1)
      call bcast_int1(writ_cn,1)
      call bcast_int1(resconmax,1)      ! number of used continous files
      call bcast_int1(resconid,ncnfile)
      call bcast_char1(rescon,ncnfile)   ! file names
      call bcast_int1(writ_mn,1)
      call bcast_int1(nkeisu,1)
      call bcast_int(ut_constituent,tcm)
      call bcast_char1(filepath_3d,1)
      call bcast_char1(datadir,1)

      call bcast_logical1(set_start)
      if( set_start) call bcast_Tdate(date_start)

#ifndef CLIMAT
      call bcast_int1(nstart_year,1)

      cal_year = nstart_year
      if(mod(cal_year,4)==0) then
         nday=366
      else
         nday=365
      endif
#else
       cal_year = 1
#ifdef NCEPSF
      nday=365
#else
      nday=360
#endif
#endif


! tentatibe
      call bcast_int1(step_3d_year,1)
      call bcast_int1(step_3d_month,1)
      call bcast_int1(step_3d_day,1)
      call bcast_int1(step_3d_hour,1)
      call bcast_int1(step_3d_min,1)
      call bcast_int1(step_3d_sec,1)

      ns_year = step_3d_year
      ns_month = step_3d_month
      ns_day = step_3d_day
      ns_hour = step_3d_hour
      ns_min = step_3d_min
      ns_sec = step_3d_sec

#ifdef BDOUT
      nwrod=nxday*4
#endif
#if defined(SGOUT) && !defined(EXP2006) && !defined(EXP2007)
      nwrsg=nxday*1
#endif
      if(ip==imaster) then
         write(*,*)' '
         write(*,*)' '
         write(*,*)'np=',np,'jm=',jm,'imaster=',ip
         write(*,*)'nfirst=',nfirst
         write(*,*)'step_3d=',step_3d,'writ_db=',writ_db,
     &        'writ_3d=',writ_3d,'nkeisu=',nkeisu,'writ_mn=',writ_mn
#ifdef BNDTIDE
     &  ,'writ_2d=',writ_2d,'writ_en=',writ_en
#endif
#ifdef CLIMAT
            write(*,*) 'cal_year=',cal_year,'nday=',nday
#endif
#ifdef BDOUT
         write(*,*)'nwrod=',nwrod
#endif
#if defined(SGOUT) && !defined(EXP2006) && !defined(EXP2007)
         write(*,*) 'nwrsg=',nwrsg
#endif
      endif
#else ! case ASSIM
      step_3d = nt*nmodf
      nfirst = 0
      writ_3d = nmodf
      necount = 0

#ifndef CLIMAT
      cal_year = 2003
      if(mod(cal_year,4)==0) then
         nday=366
      else
         nday=365
      endif
#else ! case CLIMAT
       cal_year = 1
#ifdef NCEPSF
      nday=365
#else
      nday=360
#endif
#endif !end CLIMAT
      nxday=60*60*24/nint(dtuv)
      nxyear=nday*60*60*24/nint(dtuv)
      nxmonth=nxyear/12

#endif !end ASSIM

      dtuv2=dtuv+dtuv
      dtts2=dtts+dtts
      dttr2=dttr+dttr
      ntrop=nint(dtuv/dttr)
#ifdef GL11M
      modmtn = 12
      modmtnb = 4
#endif
#ifdef EQ100M
      modmtn=12
      modmtnb=6
      writ_3d=nx10minu
      nwr1d=nx10minu
#endif
#ifdef PC68M
      modmtn = 12
!      modmtnb = 5
      modmtnb = 3
#endif
#ifdef JP44
      modmtn = 12
      modmtnb = 5
#endif
#ifdef NWNPAC
      modmtn = 12
      modmtnb = 5
#endif
#ifdef JP68M
      modmtn = 12
      modmtnb = 5
#endif
#ifdef ROKKA
      modmtn = 12
      modmtnb = 5
#endif

      if(ip.eq.imaster) then
         write(*,*)'dtuv=',dtuv,'dtts=',dtts,'dttr=',dttr
         write(*,*)'ntrop=',ntrop,'modmtn=',modmtn,'modmtnb=',modmtnb
         write(*,*)' '
      endif

!     hupp : ratio of horizontal up-current flux to total
!     vupp : ratio of vertical up-current flux to total

      hupp = 0.5d0
      vupp = 0.5d0

!     adtuv: time of one step in hour for momentum equations.
!     adtts: time of one step in hour for tracer equations.
!     adttr: time of one step in hour for vorticity equations.
      adtts=dtts/60.d0/60.d0

!cc   instrument number
!
!     topodt
!     bnddt
!     inidt
!     ndgdt
!     result
!     continu
!     restart
!
! c     topodt=50
! c     bnddt=51
! c     inidt=52
! c     ndgdt=53
! c     result=60
! c     restart=80
! c     resmean=65
! c     heatmean=66
! c     result1d=67
! c     nwrisop=68
! c     spinupf=90
! c     tracerf=71
! c     tracont=70
      topodt=11!50
      bnddt=12!51
      inidt=13!52
      ndgdt=14!53
      result=31!60
      continu=21!70
      restart=20!80
      resmean=41!65
#ifdef BNDTIDE
      recfile_2d=42 !***
      recfile_en=43 !***
#endif
      heatmean=66
      result1d=67
      nwrisop=51!68
      spinupf=63
      tracerf=71
      tracont=70
      flxmean=80
#ifdef NESTED
!   original model dataunit number
      omtopo=15
      omunt=23
#ifndef JP44
!  initial condition (pc68m continue file)
      contini=37
#endif
#endif
#ifdef ASSIM
! true data 16
! initial guess 17
! weight param 18,19
! restart 21
! continue 22
#endif
#ifdef BDOUT
      fuoutdt=81
#endif
#if defined(SGOUT) && !defined(EXP2006) && !defined(EXP2007)
      fuoutsg=83
#endif
#ifdef ROKKA
      futidemix=74
#endif
#if defined(EXP2006)&& defined(NWNPAC)
      outini=38
#endif
#ifdef MASKLVIS
      numlvis=75
#endif
c
c
ccc   parameter associated with horizontal resolution
c
c     slat  : southernmost latitude for j=1
c     dxdeg : zonal resolution in degree
c     dydeg : meridional resolution in degree
c     hdts  : horizontal diffusion coefficient in cgs unit
c     hduv  : horizontal viscosity coefficient in cgs unit
c         hduv is used for barotropic mode in case of iso-dif
c
#ifdef GL11M
      slat=-77.d0
      dxdeg=1.d0
      dydeg=1.d0

      hdts=1.d7
      hduv=1.d8

      hdts_btm=1.d7
      hdts_sfc=1.d7
#endif
#ifdef EQ100M
      slat=-0.5-2./100.
!      slat=-5./100.-2./100.
      dxdeg=1./100.
      dydeg=1./100.

      hdts=1.d5
      hduv=1.d6
#endif
#ifdef PC68M
      slat=-10.d0-2.d0/8.d0
      dxdeg=1.d0/6.d0
      dydeg=1.d0/8.d0

      hdts=1.d6
!      hduv=1.d7
      hduv=1.d6
      hdts_btm=1.d6
      hdts_sfc=1.d6
#endif
#ifdef JP44
      hdts=1.d6
      hduv=2.d6
      hdts_btm=1.d7
#endif
#ifdef ROKKA
      ! hdts=1.d6
      ! hduv=5.d6
      hdts=1.d6
      hduv=1.d6
      hdts_btm=1.d6
      hdts_sfc=1.d6
#endif
#ifdef NWNPAC
      hdts=1.d6
      hduv=1.d6
      hdts_btm=1.d6
      hdts_sfc=1.d6
#endif
#ifdef JP68M
      hdts=1.d6
      hduv=1.d6
      hdts_btm=1.d6
      hdts_sfc=1.d6
#endif

#ifdef NESTED

#ifdef CYCLIC
      write(*,*)'Do not take CYCLIC and NESTED at the same time ;('
      stop
#endif
!   slon : east-most of this model
!   slat0,slon0 : south-most and east-most of the original model
!  obc_clim : ==1  climatological open bounday condition

#ifdef JP44
!      dtbd=60.d0*75.d0

      slat0=-77.d0
      slon0=28.d0
      dxdeg0=1.d0
      dydeg0=1.d0

      obc_clim=1
c   restoring time scale

      chfht=1.d0/60.d0**2/24.d0/1.d-100
      chfhcl=chfht
      chft=1.d0/60.d0**2/24.d0/1.d-1
      chfice=1.d0/60.d0**2/24.d0/1.d-1
      chfq=1.d0/60.d0**2/24.d0/1.d0

c   interpolation parameter
c            nrng : extraporation number
c            nsm  : smoothing times
c            cay  : =0:laplace   ->inf:spline
      NRNG = 100
      NSM =  5
      CAY = 5.d0
#endif
#ifdef NWNPAC
!      dtbd=60.d0*8.d0*365.d0/360.d0

      dxdeg0=1.d0/6.d0
      dydeg0=1.d0/8.d0
!case PC68M for obc
      slat0=-10.d0-2.d0/8.d0
      slon0=105.d0
! case JP68M for obc
      slat0=-10.d0-2.d0/8.d0+dydeg0*dble(288-1)-dydeg0*3.d0
      slon0=105.d0+dxdeg0*dble(108-1)-dxdeg0*3.d0

      obc_clim=1

      chfht=1.d0/60.d0**2/24.d0/1.d-1
      chfhcl=1.d0/60.d0**2/24.d0/1.d-1
      chft=1.d0/60.d0**2/24.d0/1.d-1
      chfice=1.d0/60.d0**2/24.d0/1.d-1
      chfq=0.d0

      NRNG = 100
      NSM = 10
      CAY = 5.d0
#endif
#ifdef JP68M
      slat0=-10.d0-2.d0/8.d0
      slon0=105.d0
      dxdeg0=1.d0/6.d0
      dydeg0=1.d0/8.d0

      obc_clim=1

      chfhcl=1.d0/60.d0**2/24.d0/1.d-1
      chft=1.d0/60.d0**2/24.d0/1.d-1
      chfice=1.d0/60.d0**2/24.d0/1.d-1
      chfq=0.d0

      NRNG = 100
      NSM = 0
      CAY = 5.d0
#endif
#ifdef ROKKA
      slat0=30.5
      slon0=127.+1./6.
      dxdeg0=1.d0/18.d0
      dydeg0=1.d0/24.d0

      obc_clim=1

!     chfhcl=1.d0/60.d0**2/24.d0/1.d-1 !original
#ifdef BNDTIDE
      chfht=1.d0/60./5.
#endif
      chfhcl=chfht
      chft=1.d0/60.d0**2/24.d0/1.d-1
      chfice=1.d0/60.d0**2/24.d0/1.d-1
      chfq=0.d0

      NRNG = 100
      NSM = 3
      npresm=0
      CAY = 5.d0
#endif

      dxdeg=dxdeg0/dble(inc)
      dydeg=dydeg0/dble(jnc)

      slat=slat0+dydeg0*dble(jsm-1)-dydeg*3.d0
      slon=slon0+dxdeg0*dble(ism-1)-dxdeg*3.d0
      slonu=slon+5.d-1*dxdeg

      if(ip==imaster) then
         write(*,*) 'slat=',slat,'  slon=',slon
         write(*,*) 'img=',img,'  jmg=',jmg
      endif

      slatu0=slat0+5.d-1*dydeg0
      slonu0=slon0+5.d-1*dxdeg0
#endif

!cc   isopicnal diffusion

!     aip: isopicnal component
!     adp: diapicnal component
!     ar: ratio of aip to adp
!     agm: skew component for Gent-Mcwillams scheme
!     ags: coefficient of Green-Stone skew parameterization
c
#ifdef GL11M
      aip = 1.d7
!      aip = 1.d6
!      adp = 1.d-1
      adp = 1.d-2
      ar = adp/aip
      agm = 1.d7
!      agm = 1.d6
      ags=0.015
      drmax=1.d-11
!      sxymax = 1.d5
!      sxymax=1.d1
      sxymax=1.d-1
      hdtsmin = 1.d4
!      hdtsmin = 1.d3
#endif
#ifdef EQ100M
      aip=1.d5
      adp=1.d-2
      ar=adp/aip
      agm=1.d5
      ags=0.015
      drmax=1.d-11
      sxymax=1.d-1
      hdtsmin=1.d3
#endif
#ifdef PC68M
      aip=5.d5
!      aip=1.d6
      adp=1.d-2
      ar=adp/aip
      agm=5.d5
!      agm=1.d6
!      ags=0.015
      drmax=1.d-11
!      sxymax=1.d-1
      sxymax=1.d-2
      hdtsmin=1.d3
#endif
#ifdef JP44
      aip = 2.d6
      adp = 1.d-2
      ar = adp/aip
      agm = 2.d6
      ags=0.015
      drmax=1.d-11
      sxymax=1.d-1
      hdtsmin = 1.d3
#endif
#ifdef ROKKA
      aip=1.d5
      adp=1.d-2
      ar=adp/aip
      agm=1.d5
!      agm=0.d0
      drmax=1.d-11
      sxymax=1.d2
!      sxymax=1.d-2
      hdtsmin=1.d2
!      hdtsmin=1.d3
#endif
#ifdef NWNPAC
!      aip=5.d5
      aip=2.d5  ! case 56
      adp=1.d-2
      ar=adp/aip
!      agm=5.d5
      agm=2.d5  ! case 56
      drmax=1.d-11
      sxymax=1.d-2
      hdtsmin=1.d3
#endif
#ifdef JP68M
      aip=1.d6
      adp=1.d-2
      ar=adp/aip
      agm=1.d6
      drmax=1.d-11
      sxymax=1.d-2
      hdtsmin=1.d3
#endif
ccc   biharmonic Smagorinsky friction shceme parameter
      hduv_bh=1.d17
      csmag=2.d0
#ifdef PC68M
      csmag=2.d0
#endif
!      csmag=5.d0
#ifdef BNDTIDE
      csmag=2.d0
#endif

      if(ip.eq.imaster) then
      write(*,*)'hdts=',hdts,'hduv=',hduv
      write(*,*)'aip=',aip,'adp=',adp,'agm=',agm
      write(*,*)' '
      endif

!   definition of layer thickness

!   set the depth interval

#if defined(GL11M) || defined(JP44)
      dz(0) = 0.
      dz(1) = 2.d3
      dz(2) = 2.d3
      dz(3) = 2.d3
      dz(4) = 2.5d3
      dz(5) = 3.d3
      dz(6) = 3.d3
      dz(7) = 3.5d3
      dz(8) = 4.5d3
      dz(9) = 5.d3
      dz(10) = 5.d3
      dz(11) = 6.d3
      dz(12) = 7.5d3
      dz(13) = 9.d3
      dz(14) = 1.d4
      dz(15) = 1.d4
      dz(16) = 1.d4
      dz(17) = 1.25d4
      dz(18) = 1.5d4
      dz(19) = 1.75d4
      dz(20) = 2.d4
      dz(21) = 2.5d4
      dz(22) = 2.5d4
      dz(23) = 3.d4
      dz(24) = 3.d4
      dz(25) = 3.d4
      dz(26) = 3.d4
      dz(27) = 3.d4
      dz(28) = 4.d4
      dz(29) = 4.d4
      dz(30) = 4.d4
      dz(31) = 4.d4
      dz(32) = 4.d4
      dz(33) = 4.d4
      dz(34) = 4.d4
#endif

#ifdef EQ100M
c      data dz/ !0:225+1
c     & 0.d0, 200*50.d0, 60.d0, 70.d0, 80.d0, 90.d0, 100.d0,
c     & 120.d0, 140.d0, 160.d0, 180.d0, 200.d0, 250.d0, 300.d0, !117.50m
c     & 350.d0, 400.d0, 500.d0, 700.d0, 10.d2, 15.d2, 20.d2, !182m
c     & 30.d2, 40.d2, 50.d2, 70.d2, 100.d2, 100.d2, 0.d0/ !572m
      dz(0)=0.d0
      do k=1,200
         dz(k)=50.d0
      enddo
      dz(201)=60.d0
      dz(202)=70.d0
      dz(203)=80.d0
      dz(204)=90.d0
      dz(205)=100.d0
      dz(206)=120.d0
      dz(207)=140.d0
      dz(208)=160.d0
      dz(209)=180.d0
      dz(210)=200.d0
      dz(211)=250.d0
      dz(212)=300.d0
      dz(213)=350.d0
      dz(214)=400.d0
      dz(215)=500.d0
      dz(216)=700.d0
      dz(217)=10.d2
      dz(218)=15.d2
      dz(219)=20.d2
      dz(220)=30.d2
      dz(221)=40.d2
      dz(222)=50.d2
      dz(223)=70.d2
      dz(224)=100.d2
      dz(225)=100.d2
      dz(226)=0.d0
#endif

#if defined(PC68M) || defined(NWNPAC) || defined(ROKKA) || defined(JP68M)
      dz(0)=0.d0
      do k=1,5
         dz(k)=4.d2
      enddo
      do k=6,41
         dz(k)=5.d2
      enddo
      dz(42)=6.d2
      dz(43)=7.d2
      dz(44)=8.d2
      dz(45)=9.d2
      do k=46,52
         dz(k)=10.d2
      enddo
      do k=53,62
         dz(k)=20.d2
      enddo
      do k=63,64
         dz(k)=50.d2
      enddo
      dz(65)=100.d2
      do k=66,67
         dz(k)=200.d2
      enddo
      do k=68,69
         dz(k)=300.d2
      enddo
      do k=70,71
         dz(k)=400.d2
      enddo
      do k=72,78
         dz(k)=500.d2
      enddo
      dz(79)=0.d0
#endif
c
c
ccc   definition of levels
c
c     dzz : distance between vertically adjascent grid levels
c     dp : depth of grid levels
c     dep : depth of layer boundary
c     dzr,dzzr : reversed values of dz and dzz
c
      dzz(1) = 0.5*dz(1)
      do k = 2,km
        dzz(k) = 0.5*(dz(k)+dz(k-1))
      enddo
      dzz(km+1) = dz(km)*0.5
      dep(0)=0.d0
      dp(1)=dzz(1)
      dep(1)=0.d0
      dep(2)=dz(1)
      dzr(1)=1.d0/dz(1)
      dzzr(1)=1.d0/dzz(1)
      do k=2,km
        dp(k)=dp(k-1)+dzz(k)
        dep(k+1)=dep(k)+dz(k)
        dzr(k)=1.d0/dz(k)
        dzzr(k)=1.d0/dzz(k)
      enddo
      dzzr(km+1)=1.d0/dzz(km+1)

!cc   definition of constants

!     omega   : rotation rate of Earth
!     vdts    : coefficient for vertical diffusion (z-dependent)
!     dxddy   : area of horizontal grid cell on sphere

      pi=3.141592653589793d0
      omega=pi/43082.d0
      radian=180.d0/pi
      radius=6375.d5
      grav=981.d0

      dx=dxdeg*radius/radian
      dy=dydeg*radius/radian
      dxr=1.d0/dx
      dyr=1.d0/dy
      dx2r=0.5d0*dxr
      dy2r=0.5d0*dyr
      dx4r=0.25d0*dxr
      dy4r=0.25d0*dyr
      dxdyr=dx*dyr
      dydxr=dy*dxr
      dxsq=dx*dx
      dysq=dy*dy
      dxsqr=1.d0/dxsq
      dysqr=1.d0/dysq
      ddy=2.d0*radius*dsin(.5d0*dydeg/radian)
      tanfi4=dtan(.25d0*dydeg/radian)
      dxddy=dx*ddy
      slatu=slat+0.5d0*dydeg
      stan=1.d-15
#ifdef EQ100M
      ddy=dy
      dxddy=dx*ddy
#endif
c
ccc   these arrays are local
c
      do j=1,jm
        jlc = ls(ip_y) + j -1
#ifdef EQ100M
         cs(j)=1.d0
         sine(j)=0.d0
         tng(j)=0.d0
         csr(j)=1.d0
         cst(j)=1.d0
         cor(j)=2.d0*omega/radius
     &        *(slatu+dble(jlc-1)*dydeg)*radius/radian
#else
         cs(j)=dcos((slatu+dble(jlc-1)*dydeg)/radian)
         sine(j)=dsin((slatu+dble(jlc-1)*dydeg)/radian)
         tng(j)=sine(j)/cs(j)
         csr(j)=1.d0/cs(j)
         cst(j)=dcos((slat+dble(jlc-1)*dydeg)/radian)
         cor(j)=2.d0*omega*sine(j)
#endif
         ashf(j)=.25d0*(1.d0+tng(j)*tanfi4)*cs(j)*dxddy
         anhf(j)=.25d0*(1.d0-tng(j)*tanfi4)*cs(j)*dxddy
         areauu(j)=dxddy*cs(j)
         areaur(j)=1.d0/areauu(j)
         cstr(j)=1.d0/cst(j)

#ifdef BISMFRIC
         delx(j) = dmin1(dx*cs(j),dy)
#endif
          enddo
c
ccc   parameter for bottom stress
c
      bsn=dsin(1.d1/radian)
      bcs=dcos(1.d1/radian)
      lateq = nint(-slat/dydeg)
      do j=1,jm
         jlc = ls(ip_y) + j -1
         if(jlc.le.lateq) then
            sgn(j)=-1.d0*bsn
         else
            sgn(j)=1.d0*bsn
            endif
         enddo
c
ccc   parameters for mixed layer model
c
c     ckrmn: von Karman constant
c     rhogh0: rhoughness length scale at sea surface
c     calph: constant factor for Ri
c     vdtsmn: minimum value of vdts
c     vduvmn: minimun value of vduv
c     vdtsmx: maximum value of vdts
c     c0: contant for dissipation
c     s0: constant for viscosity
c     cpr: ratio of tracer diffusivity to viscosity
c     sigma: ratio of tke diffusivity to viscosity
c     tkemin: minimum tke
c     drfit: fitting coeficients of drho(sst)
c
      ckrmn = 0.4d0
      rough0 = 1.d2
      calph = 100.d0
      cm0  = 500.d0
      vdtsmn = 1.d-2
      vduvmn = 2.d-1
      vdtsmx = 1000.d0
      vduvmx = 1000.d0
      c0 = 0.06d0
      s0 = 0.39d0
      cpr = 0.8d0
      sigma = 1.95d0
      tkemin = 1.d-6

      dleng_e=3.d2
      c0_e=1.d0
      calph_e=0.082d0
      pr0=0.2d0
      calph_p=0.82d0
      calph_b=5.27d0

      drfit(1) =  5.90257d-05
      drfit(2) = -8.60480d-06
      drfit(3) =  1.23469d-06
      drfit(4) = -6.03051d-08
      drfit(5) =  1.33955d-09
c
#ifdef PC68M
      vdtsmx=500.d0
      vduvmx=500.d0
#endif
#ifdef ROKKA
      vdtsmx=500.d0
      vduvmx=500.d0
#endif
#ifdef NWNPAC
      vdtsmx=500.d0
      vduvmx=500.d0
#endif
#ifdef JP68M
      vdtsmx=500.d0
      vduvmx=500.d0
#endif
#ifdef JP44
      vdtsmx=500.d0
      vduvmx=500.d0
#endif
c
ccc   parameter for bulk formula
c
c     unit is mks
c     rho_a/w/i  : air/water/ice density [kg m-3]
c     c_p_a/w/i  : specific heat of air/water/ice [Ws kg-1 K-1]
c     l_w/i      : latent heat of vaporization/sublimation
c     a-d_w/i    : constant for water vaper pressure
c     eps_w      : emissivity of water
c     cboltz     : Stefan-Boltzmann constant
c     p_a        : air pressure
c     c_se       : transfer coefficient for sensible heat
c     c_la_w     : transfer coefficient for latent heat
c     albedo_w/i : albedo of water/ice
c     a1-3_frz   : constant for freezing point temperature
c     k_i        : conductivity of ice [W m-1 K-1]
c     hmin       : minimum value of ice hight [cm]
c     s_ice      : sea ice salinity [psu]
c     hdice      : horizontal diffusivity of sea ice [cm2 s-1]
c     cd_i_o     : water drag rate coefficient of sea ice
c
      rho_a=1.3d0
      rho_w=1.03d3
      rho_i=0.91d3
      c_p_a=1004.d0
      c_p_w=4.2d3
      c_p_i=3.36d5
      l_w=2.5d6
      l_i=2.834d6
      a_w=6.1121d0
      b_w=18.729d0
      c_w=257.87d0
      d_w=227.3d0
      a_i=6.1078d0
      b_i=23.036d0
      c_i=279.82d0
      d_i=33.7d0
      eps_w=0.97d0
      cboltz=5.5d-8
      p_a=1024.d0
      c_stab=0.83d-3  !ikeda
      c_unst=1.10d-3
      c_la_w=1.2d-3   !1.5e-3 ikeda
      albedo_i=0.7d0  !ikeda
      albedo_w=0.1d0  !ikeda
      albedo_melt=0.4d0
      a1_frz=-0.0575d0
      a2_frz=1.710523d-3
      a3_frz=-2.154996d-4
      k_i=2.04d0
      hmin=50.d0
      hmax=500.d0
      s_ice=3.d0
      hdice=1.d7
c      hdice=1.d8
      cd_a_i=2.d-3
c      cd_a_i=1.d-3
#ifdef CLIMAT
      cd_a_o=2.18d-3
#else
      cd_a_o=1.d-3
#endif
      cd_i_o=5.d-3
      hminu=1.d0
#ifdef GL11M
      dtuice=10.d0
      hminu=hmin
#endif
#ifdef JP44
      dtuice=10.d0
      hminu=hmin
#endif
#ifdef PC68M
      dtuice=4.d0
      hdice=1.d7
      hminu=10.d0
#endif
#ifdef ROKKA
c      dtuice=4.d0
      dtuice=2.d0
      hdice=1.d7
      hminu=10.d0
#endif
#ifdef NWNPAC
c      dtuice=4.d0
      dtuice=2.d0
      hdice=1.d7
      hminu=10.d0
#endif
#ifdef JP68M
      dtuice=4.d0
      hdice=1.d7
      hminu=10.d0
#endif
c      aicemin=1.d-10
      aicemin=1.d-11
c      dtskin=0.5d0
      dtskin=0.d0
c      t_ice_lim = 0.01d0
      t_ice_lim = 0.0001d0
      modice = 5
c      dump_i = 1.d-4
c      dump_i = 1.d-5
      dump_i = 1.d-6
      do j=1,jm
         xkai(j)=0.5+0.4*dabs(slat+dble(ls(ip_y)+j-1)*dydeg)/90.
      enddo
#ifdef PC68M
      do j=1,jm
         xkai(j)=0.5+0.4*dabs(slat+dble(ls(ip_y)+j-1)*dydeg)/90.
      enddo
#endif
c
c
ccc   read topography data and surface boundary data
c
c     ho4 : depth in cm (int*4(im,jm)) <u/v grid>
c     exn : levels of uv-box (int*4(im,jm))
c     wsx,wsy : wind stress (real*8(im,jm)) <u/v grid>
c     tini,sini : intial condition of t/s (real*8(im,jm,km))
c     tref,sref : reference value of t/s for gamma term
c     gref : strength of t/s gamma term (real*8(im,jm,km))
c
c
#ifdef GL11M
      call read_2d_i(ho4,topodt)
      call read_2d_i(exn,topodt)
#endif
#ifdef EQ100M
      do j=1,jml
      do i=1,im
         ho4(i,j)=nint(dep(km+1))
         exn(i,j)=km
         jj=ls(ip_y)+j-1
         if(jj.eq.1.or.jj.eq.2.or.
     &        jj.eq.jmg-2.or.jj.eq.jmg-1.or.jj.eq.jmg)then
            ho4(i,j)=0
            exn(i,j)=0
         endif
      enddo
      enddo
#endif
#ifdef PC68M
      call read_2d_i(ho4,topodt)
      call read_2d_i(exn,topodt)
#endif
#ifdef JP44
      call read_2d_i(ho4,topodt)
      call read_2d_i(exn,topodt)
#endif
#ifdef ROKKA
      call read_2d_i(ho4,topodt)
      call read_2d_i(exn,topodt)
#endif
#ifdef NWNPAC
      call read_2d_i(ho4,topodt)
      call read_2d_i(exn,topodt)
#endif
#ifdef JP68M
      call read_2d_i(ho4,topodt)
      call read_2d_i(exn,topodt)
#endif
#ifdef GL11M
        call read_ref(tref12,inidt)
        call read_ref(sref12,inidt)
#endif
#ifdef PC68M
        do m=1,imn
        call read_ref1(tsref1,inidt)
        do j=1,jml
        do i=1,iml
        do k=0,km+1
           tref12(i,j,k,m)=tsref1(i,j,k)
        enddo
        enddo
        enddo
        enddo
        do m = 1,imn
        call read_ref1(tsref1,inidt)
        do j=1,jml
        do i=1,iml
        do k=0,km+1
           sref12(i,j,k,m)=tsref1(i,j,k)
        enddo
        enddo
        enddo
        enddo
#endif
#ifdef JP68M
        do m=1,imn
        call read_ref1(tsref1,inidt)
        do j=1,jml
        do i=1,im
        do k=0,km+1
           tref12(i,j,k,m)=tsref1(i,j,k)
        enddo
        enddo
        enddo
        call read_ref1(tsref1,inidt)
        do j=1,jml
        do i=1,im
        do k=0,km+1
           sref12(i,j,k,m)=tsref1(i,j,k)
        enddo
        enddo
        enddo
        enddo
#endif

!   read nadging factor
#ifdef GL11M
      call read_3d_r48(gref,ndgdt)
#endif
#ifdef EQ100M
      do j=1,jml
      do k=1,km
      do i=1,im
         gref(i,j,k)=0.
      enddo
      enddo
      enddo
#endif
#ifdef PC68M
      call read_3d_r4(gref,ndgdt)

!      do j = 1,jml
!      do i = 1,iml
!        gref(i,j,1) = 1.d0/(60.d0*60.d0*24.d0*10.d0)
!      enddo
!      enddo

!      do k = 2,km
!      do j = 1,jml
!      do i = 1,iml
!        gref(i,j,k) = 1.d0/(60.d0*60.d0*24.d0*30.d0)
!      enddo
!      enddo
!      enddo

#endif
#ifdef JP44
      do j=1,jml
      do k=1,km
      do i=1,im
         gref(i,j,k)=0.
      enddo
      enddo
      enddo
#endif
#ifdef ROKKA
      do j=1,jml
      do k=1,km
      do i=1,iml
         gref(i,j,k)=1.d0/(60.d0**2*24*30.d0)
      enddo
      enddo
      enddo
#endif
#ifdef NWNPAC
#ifdef CASE58
      call read_3d_r8(gref,ndgdt)
      do j=1,jml
         grefini(j)=1.d0/(60.d0**2*24*10.d0)
      enddo
! south
      if(ip_y .eq. 0) then
         do j=1,jndg+2
            grefini(j)=0.d0
         enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
         do j=1,jndg+2
            grefini(jml-1-j)=0.d0
         enddo
      endif
#elif defined(CASE60)
      call read_3d_r48(gref,ndgdt)

      do j=1,jml
         grefini(j)=1.d0/(60.d0**2*24*10.d0)
      enddo
! south
      if(ip_y .eq. 0) then
         do j=1,jndg+2
            grefini(j)=0.d0
         enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
         do j=1,jndg+2
            grefini(jml-1-j)=0.d0
         enddo
      endif
#else
      do j=1,jml
      do k=1,km
      do i=1,iml
         gref(i,j,k)=1.d0/(60.d0**2*24*10.d0)
      enddo
      enddo
      enddo
#endif

! west
      if(ip_x .eq. 0) then
         do j=1,jml-1
         do i=1,indg+2
            gref(i,j,1:km)=0.d0
         enddo
         enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
         do j=1,jml-1
         do i=1,indg+2
            gref(iml-1-i,j,1:km)=0.d0
         enddo
         enddo
      endif
! south
      if(ip_y .eq. 0) then
         do i=1,iml-1
         do j=1,jndg+2
            gref(i,j,1:km)=0.d0
         enddo
         enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
         do i=1,iml-1
         do j=1,jndg+2
            gref(i,jml-1-j,1:km)=0.d0
         enddo
         enddo
      endif
#endif
#ifdef JP68M
      call read_3d_r4(gref,ndgdt)
#endif
c

#ifdef MASKLVIS
      call read_2d_r8(mask_lvis,numlvis)
#endif

ccc   daily forcing
c
c      read(bnddt)wsx_d
c      read(bnddt)wsy_d
c      read(bnddt)sc_wind_d
c      read(bnddt)stdev_w_d
c      read(bnddt)dpt_t2m_d
c      read(bnddt)temp_2m_d
c      read(bnddt)cloud_d
c      read(bnddt)tot_sol_d
c      read(bnddt)wflux_d
c
#ifdef GL11M
      call read_fd(wsx_d,bnddt)
      call read_fd(wsy_d,bnddt)
      call read_fd(sc_wind_d,bnddt)
      call read_fd(stdev_w_d,bnddt)
      call read_fd(dpt_t2m_d,bnddt)
      call read_fd(temp_2m_d,bnddt)
      call read_fd(cloud_d,bnddt)
      call read_fd(tot_sol_d,bnddt)
      call read_fd(wflux_d,bnddt)
#endif

#ifdef EQ100M
      if(ip.eq.imaster)then
      read(bnddt)awind_h
      read(bnddt)xwind_h
      read(bnddt)ywind_h
      read(bnddt)atemp_h
      read(bnddt)rhmdt_h
      read(bnddt)preci_h
      read(bnddt)apres_h
      read(bnddt)swrad_h
!      write(*,*)awind_h(1),xwind_h(1),ywind_h(1)
!      write(*,*)atemp_h(1),rhmdt_h(1),preci_h(1)
!      write(*,*)apres_h(1),swrad_h(1)
      endif
      call bcast_real(awind_h,nhour)
      call bcast_real(xwind_h,nhour)
      call bcast_real(ywind_h,nhour)
      call bcast_real(atemp_h,nhour)
      call bcast_real(rhmdt_h,nhour)
      call bcast_real(preci_h,nhour)
      call bcast_real(apres_h,nhour)
      call bcast_real(swrad_h,nhour)
!      if(ip.eq.0)then
!      write(*,*)ip,awind_h(1),xwind_h(1),ywind_h(1)
!      write(*,*)ip,atemp_h(1),rhmdt_h(1),preci_h(1)
!      write(*,*)ip,apres_h(1),swrad_h(1)
!      else if(ip.eq.1)then
!      write(*,*)ip,awind_h(1),xwind_h(1),ywind_h(1)
!      write(*,*)ip,atemp_h(1),rhmdt_h(1),preci_h(1)
!      write(*,*)ip,apres_h(1),swrad_h(1)
!      endif
!      stop
#endif

#ifdef PC68M
#ifndef NCEPSF
      call read_fd(wsx_d,bnddt)
      call read_fd(wsy_d,bnddt)
      call read_fd(sc_wind_d,bnddt)
      call read_fd(stdev_w_d,bnddt)
      call read_fd(dpt_t2m_d,bnddt)
      call read_fd(temp_2m_d,bnddt)
      call read_fd(cloud_d,bnddt)
      call read_fd(tot_sol_d,bnddt)
      call read_fd(wflux_d,bnddt)
#endif
#endif


#ifdef NCIO
      if( ip == imaster) then
#ifdef ASSIM
       write(*,*) 'write on : ',rsdir,workdir
       write(*,*) 'prev case : ',trim(incase),
     &  ', current case : ',trim(rscase)
#else
        write(*,*) 'result file : ',rsfile
#if defined(ROKKA) && defined(EXP2007)
        do n=1,step_3d/nxday
          ! write(*,*) 'mean file : ',rmfile(n)
        enddo
#else
        ! write(*,*) 'mean file : ',rmfile
#endif
#ifdef BNDTIDE
        write(*,*) 'tide file : ',filename_2d !***
#endif
#endif
      endif
#ifdef JP44
      call create_nc_result(result,rsfile,
     $     'Near Japan 1/4 x 1/4deg ocean model (snap shot)')
!      call create_nc_result(resmean,rmfile,
!     $     'Near Japan 1/4 x 1/4 deg ocean model (monthly mean)')
      call create_nc_flux(flxmean,flxmfile,
     $   'Surface flux of near Japan 1/4 x 1/4 deg ocean model'//
     $     ' (monthly mean)')
#endif
#ifdef JP68M
#ifdef ASSIM
      if(ip==imaster) then
        fwdfile=trim(workdir)//'/fwd_'//trim(rscase)//'.nc'
        write(*,*) 'forward file : ',fwdfile
        fcstfile=trim(rsdir)//'/fcst_'//trim(rscase)//'.nc'
        write(*,*) 'result file : ',fcstfile
        fcinifile=trim(rsdir)//'/fcinit_'//trim(rscase)//'.nc'
        write(*,*) 'fcinit file : ',fcinifile
      endif

      if(ip==imaster) then
      outfile = trim(rsdir)//'/last_'//trim(rscase)//'.dat'
      write(*,*) 'continue file : ',outfile
      open(continu,outfile,form='unformatted')

      outfile = trim(rsdir)//'/rslt_'//trim(rscase)//'.dat'
      write(*,*) 'result file (obs space) : ',outfile
      open(32,outfile,form='unformatted')

      outfile = trim(rsdir)//'/var_'//trim(rscase)//'.dat'
      write(*,*) 'control variable file : ',outfile
      open(40,outfile,form='unformatted')

      outfile = trim(rsdir)//'/grd_'//trim(rscase)//'.dat'
      write(*,*) 'gradient file : ',outfile
      open(41,outfile,form='unformatted')

      outfile = trim(rsdir)//'/dir_'//trim(rscase)//'.dat'
      write(*,*) 'direction of iteration file : ',outfile
      open(42,outfile,form='unformatted')

      outfile = trim(rsdir)//'/init_'//trim(rscase)//'.dat'
      write(*,*) 'initilal condition file : ',outfile
      open(43,outfile,form='unformatted')

      outfile = trim(rsdir)//'/it'//trim(rscase)//'.dat'
      write(*,*) 'iteration output file : ',outfile
      open(22,outfile,form='unformatted')

      if(rscase/=incase) then
        outfile = trim(rsdir)//'/it'//trim(incase)//'.dat'
        write(*,*) 'iteration input file : ',outfile
        open(24,outfile,action='READ',form='unformatted')
      endif

      endif
! without netcdf flux
      call read_fd(wsx_d,bnddt)
      call read_fd(wsy_d,bnddt)
      call read_fd(sc_wind_d,bnddt)
      call read_fd(stdev_w_d,bnddt)
      call read_fd(dpt_t2m_d,bnddt)
      call read_fd(temp_2m_d,bnddt)
      call read_fd(cloud_d,bnddt)
#ifdef NCEPSF
      call read_fd(dswrf_d,bnddt)
#else
      call read_fd(tot_sol_d,bnddt)
#endif
      call read_fd(wflux_d,bnddt)

#else ! ASSIM
      call create_nc_result(result,rsfile,
     $     'Near Japan 1/6 x 1/8 deg ocean model (snap shot)')
     ! call create_nc_result(resmean,rmfile,
     ! $     'Near Japan 1/6 x 1/8 deg ocean model (monthly mean)')
      call create_nc_flux(flxmean,flxmfile,
     $   'Surface flux of near Japan 1/6 x 1/8 deg ocean model'//
     $     ' (monthly mean)')
#endif ! end ASSIM
#endif ! end PC68M
#ifdef NWNPAC
      call create_nc_result(result,rsfile,
     $     'Near Japan 1/18 x 1/24 deg ocean model (snap shot)')
     ! call create_nc_result(resmean,rmfile,
     ! $     'Near Japan 1/18 x 1/24 deg ocean model (monthly mean)')
      call create_nc_flux(flxmean,flxmfile,
     $   'Surface flux of near Japan 1/18 x 1/24 deg ocean model'//
     $     ' (monthly mean)')
#endif
#ifdef ROKKA
      call create_nc_result(result,rsfile,
     &     'Off Rokkasho ocean model (snap shot)')
#ifdef EXP2007
      ! do n=1,step_3d/nxday
#ifdef DEBUG
      !    if(n>numrm) then
      !       if(ip==imaster) then
      !          write(*,*)'"numrm" must be larger !'
      !          stop
      !       endif
      !    endif
#endif
      !    call create_nc_result(resmean(n),rmfile(n),
      ! &        'Off Rokkasho ocean model')
      ! enddo
#else
      ! call create_nc_result(resmean,rmfile,
      ! &     'Off Rokkasho ocean model (monthly mean)')
#endif
      ! call create_nc_flux(flxmean,flxmfile,
      ! &   'Surface flux of off Rokkasho ocean model'//
      ! &     ' (monthly mean)')

#ifdef BNDTIDE
      call create_nc_2d(recfile_2d,filename_2d,    !***
     &     'Rokkasho short-period tide record ')

      call create_nc_en(recfile_en,filename_en,    !***
     &     'Rokkasho energy state record ')

#ifdef IDEBUG
      call create_nc_eqterms(recfile_eq,filename_eq,    !***
     &     'Balanced terms in momentum equations')
#endif


#endif

#endif


!------------------ input surface flux --------------------------------
#if defined(ROKKA) && (defined(EXP2006) || defined(EXP2007))

      do l=0,2

         if(ip.eq.imaster) then
            select case(l)
            case(0)
               ncstat=nf90_open(sffile0,NF90_NOWRITE,bnddt)
            case(1)
               ncstat=nf90_open(sffile1,NF90_NOWRITE,bnddt)
            case(2)
               ncstat=nf90_open(sffile2,NF90_NOWRITE,bnddt)
            end select
#ifdef DEBUG
            if(ncstat.ne.0) then
               write(*,*) nf90_strerror(ncstat)
               stop
            endif
#endif

            ncstat=nf90_inq_varid(bnddt,'nsec',fdid)
#ifdef DEBUG
            if(ncstat.ne.0) then
               write(*,*) nf90_strerror(ncstat)
               stop
            endif
#endif
            ncstat=nf90_get_var(bnddt,fdid,nsec_scf(l*nsf+1:(l+1)*nsf))
#ifdef DEBUG
            if(ncstat.ne.0) then
               write(*,*) nf90_strerror(ncstat)
               stop
            endif
#endif
#ifdef DEBUG13
            if(ip==imaster) then
               do nn=1,nsf
                  write(*,"('nsec_scf=',i)") nsec_scf(l*nsf+nn)
               enddo
            endif
#endif
         endif

         call read_nc_fd(sc_wind_d(1,1,l*nsf+1),bnddt,1,nsf,'sc_wind')
         call read_nc_fd(dpt_t2m_d(1,1,l*nsf+1),bnddt,1,nsf,'dpt_t2m')
         call read_nc_fd(temp_2m_d(1,1,l*nsf+1),bnddt,1,nsf,'temp_2m')
         call read_nc_fd(dswrf_d(1,1,l*nsf+1),bnddt,1,nsf,'dswrf')
#ifdef DLWRF
         call read_nc_fd(dlwrf_d(1,1,l*nsf+1),bnddt,1,nsf,'dlwrf')
#else
         call read_nc_fd(cloud_d(1,1,l*nsf+1),bnddt,1,nsf,'cloud')
#endif
         call read_nc_fd(wflux_d(1,1,l*nsf+1),bnddt,1,nsf,'wflux')
#ifndef EVAPO
         call read_nc_fd(evapo_d(1,1,l*nsf+1),bnddt,1,nsf,'evapo')
         do n=l*nsf+1,(l+1)*nsf
         do j=1,jml
         do i=1,iml
            wflux_d(i,j,n)=wflux_d(i,j,n)-evapo_d(i,j,n)
         enddo
         enddo
         enddo
#endif
         if(ip.eq.imaster) ncstat=nf90_close(bnddt)

! --- for wind stress

         if(ip.eq.imaster) then
         select case(l)
         case(0)
            ncstat=nf90_open(wsfile0,NF90_NOWRITE,bnddt)
         case(1)
            ncstat=nf90_open(wsfile1,NF90_NOWRITE,bnddt)
         case(2)
            ncstat=nf90_open(wsfile2,NF90_NOWRITE,bnddt)
         end select
#ifdef DEBUG
            if(ncstat.ne.0) then
               write(*,*) nf90_strerror(ncstat)
               stop
            endif
#endif

            ncstat=nf90_inq_varid(bnddt,'nsec',fdid)
#ifdef DEBUG
            if(ncstat.ne.0) then
               write(*,*) nf90_strerror(ncstat)
               stop
            endif
#endif

            ncstat=nf90_get_var(bnddt,fdid,nsec_ws(l*nws+1:(l+1)*nws))

#ifdef DEBUG
            if(ncstat.ne.0) then
               write(*,*) nf90_strerror(ncstat)
               stop
            endif
#endif
#ifdef DEBUG13

            do nn=1,nws
               write(*,"('nsec_ws=',i)") nsec_ws(l*nws+nn)
            enddo

#endif
         endif


         call read_nc_fd(wsx_d(1,1,l*nws+1),bnddt,1,nws,'wsx')
         call read_nc_fd(wsy_d(1,1,l*nws+1),bnddt,1,nws,'wsy')


         if(ip.eq.imaster) ncstat=nf90_close(bnddt)

      end do

      call bcast_int(nsec_scf,nsf*3)
      call bcast_int(nsec_ws,nws*3)
#else
!------------------- last day of previous year ---------------------
#ifndef ROKKA
#ifndef CLIMAT
      if(ip.eq.imaster) then
         ncstat=nf90_open(sffile0,NF90_NOWRITE,bnddt)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
#if defined(EXP2006) || defined(EXP2007)
         ncstat=nf90_inq_dimid(bnddt,'time',dimid_day)
#else
         ncstat=nf90_inq_dimid(bnddt,'day',dimid_day)
#endif
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
         ncstat=nf90_inquire_dimension(bnddt,dimid_day,len=nsfr)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
      endif
      call read_nc_fd(wsx_d(1,1,0),bnddt,nsfr,nsfr,'wsx')
      call read_nc_fd(wsy_d(1,1,0),bnddt,nsfr,nsfr,'wsy')
      call read_nc_fd(sc_wind_d(1,1,0),bnddt,nsfr,nsfr,'sc_wind')
      call read_nc_fd(dpt_t2m_d(1,1,0),bnddt,nsfr,nsfr,'dpt_t2m')
      call read_nc_fd(temp_2m_d(1,1,0),bnddt,nsfr,nsfr,'temp_2m')
#ifdef NCEPSF
      call read_nc_fd(dswrf_d(1,1,0),bnddt,nsfr,nsfr,'dswrf')
#ifdef DLWRF
      call read_nc_fd(dlwrf_d(1,1,0),bnddt,nsfr,nsfr,'dlwrf')
#else
      call read_nc_fd(cloud_d(1,1,0),bnddt,nsfr,nsfr,'cloud')
#endif

#else ! not NCEPSF
      call read_nc_fd(tot_sol_d(1,1,0),bnddt,nsfr,nsfr,'tot_sol')
      call read_nc_fd(cloud_d(1,1,0),bnddt,nsfr,nsfr,'cloud')
#endif
      call read_nc_fd(wflux_d(1,1,0),bnddt,nsfr,nsfr,'wflux')
#ifndef EVAPO
      call read_nc_fd(evapo_d(1,1,0),bnddt,nsfr,nsfr,'evapo')
      do j=1,jml
      do i=1,iml
         wflux_d(i,j,0)=wflux_d(i,j,0)-evapo_d(i,j,0)
      enddo
      enddo
#endif
      if(ip.eq.imaster) ncstat=nf90_close(bnddt)
#endif ! end CLIMAT
#endif ! end ROKKA

!------------------ 1-year flux ----------------------------------

#ifdef ROKKA
      nsfr=nday*4
#else
      nsfr=nday
#endif

      if(ip.eq.imaster) then
         ncstat=nf90_open(sffile1,NF90_NOWRITE,bnddt)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) sffile1
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
         write(*,*) 'sflux1=',sffile1
         write(*,*) 'nday=',nday,'  nsfr=',nsfr
#endif
      endif
      call read_nc_fd(wsx_d(:,:,1:nsfr),bnddt,1,nsfr,'wsx')
      call read_nc_fd(wsy_d(:,:,1:nsfr),bnddt,1,nsfr,'wsy')
      call read_nc_fd(sc_wind_d(:,:,1:nsfr),bnddt,1,nsfr,'sc_wind')
#ifdef CLIMAT
      call read_nc_fd(stdev_w_d(:,:,1:nsfr),bnddt,1,nsfr,'stdev_w')
#else
      stdev_w_d(:,:,:)=0.d0
#endif
      call read_nc_fd(dpt_t2m_d(:,:,1:nsfr),bnddt,1,nsfr,'dpt_t2m')
      call read_nc_fd(temp_2m_d(:,:,1:nsfr),bnddt,1,nsfr,'temp_2m')
#ifdef NCEPSF
      call read_nc_fd(dswrf_d(:,:,1:nsfr),bnddt,1,nsfr,'dswrf')
#ifdef DLWRF
      call read_nc_fd(dlwrf_d(:,:,1:nsfr),bnddt,1,nsfr,'dlwrf')
#else
      call read_nc_fd(cloud_d(:,:,1:nsfr),bnddt,1,nsfr,'cloud')
#endif
#else ! not NCEPSF
      call read_nc_fd(tot_sol_d(:,:,1:nsfr),bnddt,1,nsfr,'tot_sol')
      call read_nc_fd(cloud_d(:,:,1:nsfr),bnddt,1,nsfr,'cloud')
#endif !end NCEPSF
      call read_nc_fd(wflux_d(:,:,1:nsfr),bnddt,1,nsfr,'wflux')
#ifndef EVAPO
      call read_nc_fd(evapo_d(:,:,1:nsfr),bnddt,1,nsfr,'evapo')
      do n=1,nsfr
      do j=1,jml
      do i=1,iml
         wflux_d(i,j,n)=wflux_d(i,j,n)-evapo_d(i,j,n)
      enddo
      enddo
      enddo
#endif
      ncstat=nf90_close(bnddt)

!-------------- first day of next year -------------------------------
#ifndef CLIMAT
      if(ip.eq.imaster) then
         ncstat=nf90_open(sffile2,NF90_NOWRITE,bnddt)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif ! end DEBUG
      endif
      call read_nc_fd(wsx_d(1,1,nsfr+1),bnddt,1,1,'wsx')
      call read_nc_fd(wsy_d(1,1,nsfr+1),bnddt,1,1,'wsy')
      call read_nc_fd(sc_wind_d(1,1,nsfr+1),bnddt,1,1,'sc_wind')
      call read_nc_fd(dpt_t2m_d(1,1,nsfr+1),bnddt,1,1,'dpt_t2m')
      call read_nc_fd(temp_2m_d(1,1,nsfr+1),bnddt,1,1,'temp_2m')
#ifdef NCEPSF
      call read_nc_fd(dswrf_d(1,1,nsfr+1),bnddt,1,1,'dswrf')
#ifdef DLWRF
      call read_nc_fd(dlwrf_d(1,1,nsfr+1),bnddt,1,1,'dlwrf')
#else
      call read_nc_fd(cloud_d(1,1,nsfr+1),bnddt,1,1,'cloud')
#endif
#else ! not NCEPSF
      call read_nc_fd(tot_sol_d(1,1,nsfr+1),bnddt,1,1,'tot_sol')
      call read_nc_fd(cloud_d(1,1,nsfr+1),bnddt,1,1,'cloud')
#endif !end NCEPSF
      call read_nc_fd(wflux_d(1,1,nsfr+1),bnddt,1,1,'wflux')
#ifndef EVAPO
      call read_nc_fd(evapo_d(1,1,nsfr+1),bnddt,1,1,'evapo')
      do j=1,jml
      do i=1,iml
         wflux_d(i,j,nsfr+1)=wflux_d(i,j,nsfr+1)-evapo_d(i,j,nsfr+1)
      enddo
      enddo
#endif
      if(ip.eq.imaster) ncstat=nf90_close(bnddt)
#endif ! end not CLIMAT
#endif
#else
      call read_fd(wsx_d,bnddt)
      call read_fd(wsy_d,bnddt)
      call read_fd(sc_wind_d,bnddt)
      call read_fd(stdev_w_d,bnddt)
      call read_fd(dpt_t2m_d,bnddt)
      call read_fd(temp_2m_d,bnddt)
      call read_fd(cloud_d,bnddt)
#ifdef NCEPSF
      call read_fd(dswrf_d,bnddt)
#else
      call read_fd(tot_sol_d,bnddt)
#endif
      call read_fd(wflux_d,bnddt)
#endif


!   setup for adoptive coodinate

      depadp = dep(kadp+1)
      if(ip.eq.imaster) then
         write(*,*) 'level and depth of seam ',kadp,depadp
         write(*,*)' '
      endif
c
#ifdef DEBUG
      do j = 1,jml
      do i = 1,iml
      if(ho4(i,j).lt.depadp.and.ho4(i,j).gt.0) then
         write(*,*) ' topodt error (ho4<depadp): '
     $        ,ip,i,j,ho4(i,j),exn(i,j)
         ho4(i,j) = int(depadp)
         exn(i,j) = kadp
      endif
      enddo
      enddo
#endif
c
      do j=1,jml
      do i=1,iml
         hrr(i,j)=dble(ho4(i,j))
         if(exn(i,j).eq.0) then
            dzub(i,j)=0.d0
         else
            dzub(i,j)=hrr(i,j)-dep(exn(i,j))
            idep=dzub(i,j)+.0001d0
            idpst=dz(exn(i,j))*.1d0+.0001d0
            if(idep.lt.idpst) then
               dzub(i,j)=idpst
               ho4(i,j)=ho4(i,j)-idep+idpst
               hrr(i,j)=dble(ho4(i,j))
            endif
         endif
c
#ifdef DEBUG
         if(exn(i,j).eq.0)then
            if(hrr(i,j).ne.0.)
     &           write(*,*)' topodt error: ',
     &           lw(ip_x)+i-1,ls(ip_y)+j-1,hrr(i,j),exn(i,j)
         else
            if(hrr(i,j).lt.dp(exn(i,j)).or.
     &           hrr(i,j).gt.dep(exn(i,j)+1)) then
                write(*,*)' topodt error: ',
     &           lw(ip_x)+i-1,ls(ip_y)+j-1,hrr(i,j),exn(i,j)
                write(*,*) dp(exn(i,j)),hrr(i,j),dep(exn(i,j)+1)
            endif
         endif
#endif
      enddo
      enddo

      do i=1,im
      do j=1,jm
         texn(i,j)=0
         dztb(i,j)=0.d0
      enddo
      enddo
      do j=2,jml
      do i=2,iml
         texn(i,j)=max0(exn(i,j),exn(i-1,j),exn(i,j-1),exn(i-1,j-1))
         if(texn(i,j).eq.0) then
            dztb(i,j)=0.d0
         else
            dztb(i,j)=
     &           dmax1( hrr(i,j),hrr(i-1,j),hrr(i,j-1),hrr(i-1,j-1) )
     &           -dep(texn(i,j))
         endif
      enddo
      enddo

      call exch_2di_s1(texn,1)
      call exch_2d_s1(dztb,2)
      call wait_2di_s1(texn,1)
      call wait_2d_s1(dztb,2)

      call exch_2di_w1(texn,3)
      call exch_2d_w1(dztb,4)

      do j=1,jm
      do k=0,km+1
      do i=1,im
         ex(i,j,k)=0.d0
         tex(i,j,k)=0.d0
         dzu(i,j,k)=0.d0
         dzua(i,j,k)=0.d0
         dzt(i,j,k)=0.d0
         areat(i,j,k)=0.d0
      enddo
      enddo
      enddo

      call wait_2di_w1(texn,1)
      call wait_2d_w1(dztb,2)

      do k=1,km
      do j=1,jml
      do i=1,iml
         if(k.le.exn(i,j)) ex(i,j,k)=1.d0
         if(k.le.texn(i,j)) tex(i,j,k)=1.d0
         dzu(i,j,k)=dz(k)*ex(i,j,k)
         dzua(i,j,k)=dz(k)*ex(i,j,k)
         if(exn(i,j).eq.k) dzu(i,j,k)=dzub(i,j)
         if(exn(i,j).eq.k) dzua(i,j,k)=dzub(i,j)
         dzt(i,j,k) = dz(k)*tex(i,j,k)
         if(texn(i,j).eq.k) dzt(i,j,k)=dztb(i,j)
      enddo
      enddo
      enddo

      do k = 1,km
      do j = 2,jml
      do i = 2,iml

      volt(i,j,k)=tex(i,j,k)*
     &        ( (dzu(i-1,j-1,k)+dzu(i,j-1,k))*anhf(j-1)
     &        +(dzu(i-1,j,k)+dzu(i,j,k))*ashf(j) )
      voltr(i,j,k)=tex(i,j,k)/( volt(i,j,k)+(1.d0-tex(i,j,k)) )
      rar(i,j,k)=tex(i,j,k)/
     &     ( ex(i-1,j-1,k)+ex(i-1,j,k)+ex(i,j-1,k)+ex(i,j,k)
     &     +(1.d0-tex(i,j,k)) )
      volur(i,j,k)=ex(i,j,k)/
     &     ( areauu(j)*dzu(i,j,k)+(1.d0-ex(i,j,k)) )
      areat(i,j,k)=tex(i,j,k)*
     &     ( (ex(i-1,j-1,k)+ex(i,j-1,k))*anhf(j-1)
     &     +(ex(i-1,j,k)+ex(i,j,k))*ashf(j) )

!   parameter related to momentam advection

      cxn(i,j,k)=ex(i,j,k)*ex(i-1,j,k)
     $     *(ex(i,j-1,k)*ex(i-1,j-1,k)-ex(i,j-1,k)-ex(i-1,j-1,k)+3.d0)
      cxs(i,j,k)=ex(i,j-1,k)*ex(i-1,j-1,k)
     $     *(ex(i,j,k)*ex(i-1,j,k)-ex(i,j,k)-ex(i-1,j,k)+3.d0)
      cye(i,j,k)=ex(i,j,k)*ex(i,j-1,k)
     $     *(ex(i-1,j,k)*ex(i-1,j-1,k)-ex(i-1,j,k)-ex(i-1,j-1,k)+3.d0)
      cyw(i,j,k)=ex(i-1,j,k)*ex(i-1,j-1,k)
     $     *(ex(i,j,k)*ex(i,j-1,k)-ex(i,j,k)-ex(i,j-1,k)+3.d0)
      cne(i,j,k)=ex(i,j,k)*ex(i-1,j-1,k)
     $     *(3.d0-ex(i,j-1,k)-ex(i-1,j,k))
      cse(i,j,k)=ex(i-1,j,k)*ex(i,j-1,k)
     $     *(3.d0-ex(i-1,j-1,k)-ex(i,j,k))

!   index for viscosity

      if(dzu(i,j,k).ge.dzu(i-1,j,k)) then
         dzumin(i,j,k) = dzu(i-1,j,k)
         xind(i,j,k) = 1.
      else
         dzumin(i,j,k) = dzu(i,j,k)
         xind(i,j,k) = 0.
      endif
      if(dzu(i,j,k).ge.dzu(i,j-1,k)) then
         dzvmin(i,j,k) = dzu(i,j-1,k)
         yind(i,j,k) = 1.
      else
         dzvmin(i,j,k) = dzu(i,j,k)
         yind(i,j,k) = 0.
      endif

      enddo
      enddo
      enddo

      do j=2,jml
      do i=2,iml
        if(hrr(i,j).gt.hrr(i-1,j)) then
          hrumin(i,j)=hrr(i-1,j)
          aindm1(i,j)=1.
          dhru(i,j)=hrr(i,j)-hrr(i-1,j)
        else
          hrumin(i,j)=hrr(i,j)
          aindm1(i,j)=0.
          dhru(i,j)=hrr(i-1,j)-hrr(i,j)
        endif
        if(hrr(i,j).gt.hrr(i,j-1)) then
          hrvmin(i,j)=hrr(i,j-1)
          aindm2(i,j)=1.
          dhrv(i,j)=hrr(i,j)-hrr(i,j-1)
        else
          hrvmin(i,j)=hrr(i,j)
          aindm2(i,j)=0.
          dhrv(i,j)=hrr(i,j-1)-hrr(i,j)
        endif
      enddo
      enddo

      ddmna=0.d0
      ttvol=0.d0
      do k=0,km+1
         dmn(k)=0.d0
         tvol(k)=0.d0
      enddo

      do k=1,km
      do j = 3,jml-2
      do i = 3,iml-2
         tvol(k)=tvol(k)+volt(i,j,k)
      enddo
      enddo
      enddo
      call sum_all(tvol,km+2)

      do k = 1,km
         ttvol=ttvol+tvol(k)
      enddo

      if(ip.eq.imaster) then
         write(*,*)'ttvol= ',ttvol
         write(*,*)'tvol:'
         write(*,600) tvol
         write(*,*)' '
      endif

  600 format(3p10es15.5)

#ifdef ASSIM

       dum_ar(1)  = 0.
       dum_ar(2)  = 0.

       do j = 3,jml-2
       do k = 1,km
       do i = 3,iml-2
         dum_ar(1) = dum_ar(1)+volt(i,j,k)
         dum_ar(2) = dum_ar(2)+ex(i,j,k)
       enddo
       enddo
       enddo

       call sum_all(dum_ar,2)
        avevol = dum_ar(1)/dum_ar(2)
!       avevol = 1.d0
!       write(*,*) 'volume',ttvol,avevol,sumex

! make global topography
      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = ex(i,j,k)
      enddo
      enddo
      enddo

      ierr=0
      call mpi_allgather(sbuf,(im-4)*(jm-4)*km,mpi_double_precision,
     &     rbuf,(im-4)*(jm-4)*km,mpi_double_precision,comm,ierr)

      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(k-1)+(im-4)*km*(j-3)
     &        +(im-4)*(jm-4)*km*(m-1)
          ii2 = i-1 + lw(ipxx)
          jj2 = j-1 + ls(ipyy)
          if(ii2 .le. img .and. jj2 .le. jmg) then
          exg(ii2,jj2,k) = rbuf(n)
          endif
        enddo
        enddo
        enddo
      enddo

      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = tex(i,j,k)
      enddo
      enddo
      enddo

      ierr=0
      call mpi_allgather(sbuf,(im-4)*(jm-4)*km,mpi_double_precision,
     &     rbuf,(im-4)*(jm-4)*km,mpi_double_precision,comm,ierr)

      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(k-1)+(im-4)*km*(j-3)
     &        +(im-4)*(jm-4)*km*(m-1)
          ii2 = i-1 + lw(ipxx)
          jj2 = j-1 + ls(ipyy)
          if(ii2 .le. img .and. jj2 .le. jmg) then
          texg(ii2,jj2,k) = rbuf(n)
          endif
        enddo
        enddo
        enddo
      enddo
#endif


#ifdef GL11M
      do j=1,jml
      do k=0,km+1
      do i=1,iml
         if(k.eq.1)then
            gref(i,j,1)=tex(i,j,1)/(60.d0*60.d0*24.d0*60.d0)
         else
            gref(i,j,k)=gref(i,j,k)*tex(i,j,k)
         endif
      enddo
      enddo
      enddo
c
ccc   freeze point check for ref temperature
ccc   make initial condition for start from rest
c
      do j = 1,jml
      do i = 1,iml
c
      do k = 1,5
      do n=0,13
      if(tex(i,j,k).eq.1.d0)then
         t_frz=sref12(i,j,k,n)*(a1_frz
     &        +a2_frz*dsqrt(sref12(i,j,k,n))
     &        +a3_frz*sref12(i,j,k,n))
         if(tref12(i,j,k,n).lt.t_frz)tref12(i,j,k,n)=t_frz
      endif
      enddo
      enddo
      enddo
      enddo

#endif
c
#ifdef EQ100M
      if(ip.eq.imaster)then
      read(inidt)tref_h
      read(inidt)sref_h
      read(inidt)uref_h
      read(inidt)vref_h
      endif
      call bcast_real(tref_h,nhour*km)
      call bcast_real(sref_h,nhour*km)
      call bcast_real(uref_h,nhour)
      call bcast_real(vref_h,nhour)
c
      do j=1,jml
      do k=1,km
      do i=1,im
c         tini(i,j,k)=tref_h(1,k)*tex(i,j,k)
c         sini(i,j,k)=sref_h(1,k)*tex(i,j,k)
c
         tini(i,j,k)=0.
         sini(i,j,k)=0.
         do nhp=1,24*4+1
            tini(i,j,k)=tini(i,j,k)+tref_h(nhp,k)
            sini(i,j,k)=sini(i,j,k)+sref_h(nhp,k)
         enddo
         tini(i,j,k)=tini(i,j,k)/dble(24*4+1)*tex(i,j,k)
         sini(i,j,k)=sini(i,j,k)/dble(24*4+1)*tex(i,j,k)
      enddo
      enddo
      enddo
c
      do j=1,jml
      do i=1,im
         do k=2,km
            gref2(i,j,k)=0.
         enddo
         xlat=slatu+dble(ls(ip_y)+j-1)*dydeg
c         gref2(i,j,1)=(0.03-dabs(xlat))/0.03/(60.*60.*2.)*ex(i,j,1)
         gref2(i,j,1)=(0.4-dabs(xlat))/0.4/(60.*60.*2.)*ex(i,j,1)
         if(gref2(i,j,1).lt.0.)gref2(i,j,1)=0.
      enddo
      enddo
#endif
c
#ifdef PC68M
      do j=1,jm
         grefini(j)=1.d0/(60.d0*60.d0*24.d0*10.d0)
!     &        *dsin((slat+dble(ls(ip_y)+j-1)*dydeg)/radian)
!         grefini(j)=dmax1(grefini(j),
!     &        1.d0/(60.d0*60.d0*24.d0*30.d0*6.d0))
!         if(grefini(j).lt.1.d0/(60.d0*60.d0*24.d0*10.d0))
!     &        grefini(j)=1.d0/(60.d0*60.d0*24.d0*10.d0)
      enddo
#endif

!   tidal mixing

      cvdtide=0.d0
!      cvdtide=vdtsmx

      do k=1,km
      do j=1,jm
      do i=1,im
         vdtide(i,j,k)=0.d0
      enddo
      enddo
      enddo
#ifdef XTIDE
! parameterize tidal mixing for pc68 model
!      count1 = 0
!      count2 = 0
      do k = 1,km
      do j = 1,jml
      do i = 1,iml
!all region
	    if(hrr(i,j) <= 50.d2) then
           vdtide(i,j,k) = 500.d0
        elseif(hrr(i,j) <=100.d2) then
           vdtide(i,j,k) = 50.d0
        elseif(hrr(i,j) <= 150.d2) then
           vdtide(i,j,k) = 10.0
        elseif(hrr(i,j) <= 200.d2) then
           vdtide(i,j,k) = 1.d0
		endif
!        if(ls(ip_y)+j-1.ge.52*8+3) then
!          if(hrr(i,j) .le. 2.d4) then
!            vdtide(i,j,k) = 100.d0
!            count1 = count1 + 1
!          elseif(hrr(i,j) .le. 5.d4) then
!            vdtide(i,j,k) = 50.d0
!          elseif(hrr(i,j) .le. 1.d5) then
!            vdtide(i,j,k) = 10.d0
!            count2 = count2 + 1
!          endif
!        endif
      enddo
      enddo
      enddo
!      write(*,*) 'tidal mixing: ',ip,' strong: ',count1,' week: ',count2
#endif
#ifdef ROKKA
      call read_3d_r8(vdtide,futidemix)
!      do k = 1,km
!      do j = 1,jm
!      do i = 1,im
!         if(hrr(i,j) .le. 200.d2) then
!            vdtide(i,j,k) = 10.d0
!         endif
!         if(dydeg*(ls(ip_y)+j-1)+slat.ge.42.) then
!            if(hrr(i,j) .le. 2.d4) then
!               vdtide(i,j,k) = 100.d0
!            elseif(hrr(i,j) .le. 5.d4) then
!               vdtide(i,j,k) = 50.d0
!            elseif(hrr(i,j) .le. 1.d5) then
!               vdtide(i,j,k) = 10.d0
!            endif
!         endif
!      enddo
!      enddo
!      enddo

#endif

!   double diffusion parameter
      expt=0.23d0
      exps=0.76d0
      vddsmx=10.d0  !1.d0
      dratioc=1.6
      vddtmn=0.015

!   tsujino sheme parameter
      do k=1,km
         vdtsuji(k)=
     &     (dsin((dep(k)-2200.d2)/4400.d2*pi)+1.d0)/2.*2.9d0+0.1d0
         if(dep(k).gt.4400.d2)vdtsuji(k)=3.d0
      enddo
#ifdef BNDTIDE
      call loadtidepara
      call caldump
#endif

#ifdef BNDTIDE_SAVENAO
      call savenao(filepath_3d)
#endif

      return
      end subroutine set