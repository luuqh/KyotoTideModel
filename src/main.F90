!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                   c
!  Ocean Genral Circulation Model Depeloped in Kyoto University     c
!
!   $Id: main.F90 13 2008-12-12 10:53:49Z ishikawa $
!                                                                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use param
      use mod_mpi
      use mod_time
#ifdef NCIO
      use netcdf
#endif

      implicit none
#include "common.h"

      real*8 starttime,nowtime,lasttime
      integer nnfst,nnlst,iyear,imin,ihour,iday,month
      integer n_month,n_day,n_hour,n_min,n_sec,nsec_end
      real*8 aday,ayear
#ifdef NCIO
      integer ncstat,n,i,j
#endif
#ifdef BNDTIDE
      character(len=80) :: con_name
      integer:: conid,con_n
#endif

      call my_mpi_init()

!    setting parameters
      call set

!    setting initial value
      call init

!   calculation of horizontal-mean pressure
      call calpp

      if(ip.eq.imaster) then
         write(*,*)'ddmna= ',ddmna,' ddmnar= ',ddmnar
         write(*,*)'dmn:'
         write(*,600)dmn
         write(*,*)' '
  600 format(3p10f10.5)
      endif

      if(.not. nfirst) then
         nnfst=last+1
      else
         nnfst=1
#ifdef PC68M
! start from July 1
         n_month=7
         n_day = 1
         n_hour=0
         n_min =0
         n_sec =0
         call day2nsec(cal_year,nsec,n_month,n_day,
     &    n_hour,n_min,n_sec)
         nnfst = nsec/nint(dtuv)
       if(ip==imaster) then
         write(*,*) 'start nsec : ',nsec
       endif
#ifdef DEBUG
       call calender(cal_year,nsec,n_month,n_day,n_hour,n_min,n_sec)
       if(ip==imaster) then
         write(*,*) ' date: ',n_month,n_day,n_hour,n_min,n_sec
       endif
#endif !end DEBUG
#endif
      endif

      if(nday == 360) then
        step_3d = nxyear*ns_year+nxmonth*ns_month+nxday*ns_day
     &    + nxmin*ns_min + nxsec*ns_sec
      else
        call calender(cal_year,nsec,n_month,n_day,n_hour,n_min,n_sec)
        if(ip==imaster) then
          write(*,*)'start from',nsec,n_month,n_day
        endif
        n_month=n_month+ns_month
        n_day = n_day+ns_day
        n_hour= n_hour+ns_hour
        n_min= n_min+ns_min
        n_sec= n_sec+ns_sec
         call day2nsec(cal_year,nsec_end,n_month,n_day,
     &    n_hour,n_min,n_sec)
        step_3d = (nsec_end-nsec)/nint(dtuv)
        if(ip==imaster) then
          write(*,*)'step_3d=',step_3d,nsec_end,n_month,n_day
        endif
      endif


#ifdef IDEBUG
      ! step_3d = 1
#endif	  
	  
      nnlst=nnfst+step_3d-1
      nnmats=0

#ifdef NO_MPI
      call gettod(starttime)
#else
      starttime = mpi_wtime()
#endif
      nowtime = starttime


      date_now=date_start%tot_sec
#ifdef BNDTIDE
      u(:,:,:) = 0.
      v(:,:,:) = 0.
      ub(:,:,:) = 0.
      vb(:,:,:) = 0.
      hcl(:,:) = 0.
      hclb(:,:) = 0.
      um(:,:) = 0.
      umb(:,:) = 0.
      vm(:,:) = 0.
      vmb(:,:) = 0.
      wl(:,:,:) = 0.
#endif

!    start of time step

      if(ip.eq.imaster) write(*,*)'time step from/to (nfirst): ',
     &     nnfst,nnlst,'( ',nfirst,')'
      if(ip.eq.imaster) write(*,*)' '

      do 1000 nkai=nnfst,nnlst

      nnmats=nnmats+1
      nergy=0
      write_3d=0
      nkeis=0
      write_1d=0

!***
#ifdef BNDTIDE
      write_2d=0
      write_en=0
      write_cn=0
      con_n = 0

      if(mod(nnmats,writ_2d)==0)then
      if(step_2d.gt.0)then
      if((nnmats.ge.stat_2d).and.(nnmats.le.(stat_2d+step_2d)))then
        write_2d=1
      endif
      endif
      endif
      if(mod(nnmats,writ_en)==0)then
      if(step_en.gt.0)then
      if((nnmats.ge.stat_en).and.(nnmats.le.(stat_en+step_en)))then
        write_en=1
      endif
      endif
      endif
      if(writ_cn.gt.0)then
      if(mod(nnmats-stat_cn,writ_cn)==0)then
      if(step_cn.gt.0)then
      if((nnmats.ge.step_cn).and.(nnmats.le.(stat_cn+step_cn)))then
        con_n = int((nnmats-stat_cn)/writ_cn)
        ! write_cn=1
      endif
      endif
      endif
      endif


#endif

      if(mod(nnmats,writ_db)==0.or.nkai==nnlst) nergy=1

      if((nnmats.ge.stat_3d).and.(nnmats.le.(stat_3d+step_3d)))then
         if(mod(nnmats,writ_3d)==0) write_3d=1
      endif

      if(mod(nkai,nkeisu).eq.0) nkeis=1


#ifdef DEBUG
        aday=ahour/24.d0
        ayear=aday/360.d0
        iyear=ayear+.0000001d0
        imin=(ahour-iyear*8640.d0)*60.d0+.01d0
        ihour=imin/60
        imin=imin-ihour*60
        iday=ihour/24
        ihour=ihour-iday*24
        month=iday/30
        iday=iday-month*30

#ifndef NO_MPI
        lasttime = nowtime
        nowtime = mpi_wtime()
#endif
#endif

#ifdef NESTED
!      calc side bounday data number
      call sbnum
#endif
!   interpolate reference from monthly climatology
      call ref_month

!   set for matsuno scheme

      if(mod(nnmats,modmtn).eq.1) then
         c2dtuv=dtuv
         c2dtts=dtts
         c2dttr=dttr
         matsno=1
         matsn2=0
      else
         c2dtuv=dtuv2
         c2dtts=dtts2
         c2dttr=dttr2
         matsno=0
         matsn2=0
      endif


 1100 continue

!cc   make surface forcing

      call surf_bulk
!
!    start calculation in each sub-region

!    calculation of tracer and baroclinic part of velocity

      call tracli

!    calculation of barotropic part
      call brtrop

!    matsuno scheme
      if(matsno.eq.1) then
         matsno=0
         matsn2=1
         go to 1100
      endif
      ahour=ahour+adtts
      date_now=date_now + diff_dtuv

!#ifdef NESTED
      nsec=nsec+nint(dtts)
      if(nsec>nday*86400) call new_year
!#endif
!
ccc   write data on disks

!***
#ifdef BNDTIDE
      if(write_2d.eq.1.and.matsno.eq.0) then
        call write2d(recfile_2d)

#ifdef IDEBUG
      do j = 1,jml
      do i = 1,iml
      eq_ht(i,j)=hta(i,j)
      eq_uv(i,j)=hta(i+1,j+1)+hta(i,j+1)-hta(i+1,j)-hta(i,j)
      eq_id(i,j)=100.
      enddo
      enddo
      call writeq(recfile_eq)
#endif	  

      endif
      if(write_en.eq.1.and.matsno.eq.0) then
          call writeen
      endif
      if(write_cn.eq.1.and.matsno.eq.0) then
         if(con_n.gt.0)then
            con_name = trim(rescon(con_n))
            conid = resconid(con_n)
            call writecn(conid,con_name)
         endif
      endif
#endif

      if(write_3d.eq.1.and.matsno.eq.0) then
	  
          call write

		  
#ifdef NWNPAC
          call wrcontin
#endif
#ifdef PC68M
          call wrcontin
#endif
      endif

!  diagnosis of model state
      if(nergy.eq.1.and.matsno.eq.0) call mmchk

!cc   calculation of horizontal-mean pressure and coefficients in the
!cc   equation of state
      if(nkeis.eq.1.and.matsno.eq.0) call calpp

!cc   write mean state
!     if(writ_mn/=0 .and. matsno==0) call wrmean          !!in write.F
!
#ifdef EQ100M
!cc   write 1d data
      if(write_1d.ne.0.and.matsno.eq.0)
     &     call write1d         !!in write.F
#endif
#ifdef BDOUT
!c   output data for nesting model
      if(nwrod.ne.0 .and. matsno.eq.0) call odwrite
#endif

#ifdef PC68M
!   spin-up diagnosation
      if(mod(nkai,nxday).eq.0.and.matsno.eq.0) call spinup               !!in mmchk.F
#endif

 1000 continue

 1200 continue

      call wrcontin             !!in write.F

#ifdef NCIO
!!   netcdf file close
      if(ip.eq.imaster) then
          ncstat=nf90_close(result)

#ifdef BNDTIDE

          ncstat=nf90_close(recfile_2d)
          ncstat=nf90_close(recfile_en)
#endif

#if defined(ROKKA) && defined(EXP2007)
          ! do n=1,numrm
             ! ncstat=nf90_close(resmean(n))
          ! enddo
#else
          ! ncstat=nf90_close(resmean)
#endif
          ! ncstat=nf90_close(flxmean)
      endif
#endif

      call my_mpi_exit()

      stop
      END