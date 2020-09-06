! subroutines of utils for calender days

!  $Id: calender.F90 10 2008-11-08 10:07:39Z ishikawa $

      subroutine calender(cal_year,nsec,cal_month,cal_day,cal_hour,
     $     cal_min,cal_sec)

! nsec to calender day

      use param
      use mod_mpi

      implicit none
      integer,intent(in) :: nsec,cal_year
      integer,intent(out) :: cal_month,cal_day,cal_hour,cal_min,cal_sec
      integer :: mday(12),msum,m

      mday(1)=31
      if(mod(cal_year,4)==0) then
         mday(2)=29
      else
         mday(2)=28
      endif
      mday(3)=31
      mday(4)=30
      mday(5)=31
      mday(6)=30
      mday(7:8)=31
      mday(9)=30
      mday(10)=31
      mday(11)=30
      mday(12)=31

      cal_sec=mod(nsec,60)
      cal_min=mod(nsec/60,60)
      cal_hour=mod(nsec/3600,24)
      cal_day=nsec/86400+1
      msum=1
      do m=1,12
         if(cal_day>=msum .and. cal_day<msum+mday(m))then
            cal_month=m
            exit
         endif
         msum=msum+mday(m)
      enddo
      cal_day=cal_day-msum+1


      return
      end subroutine calender

! ----------------------------------------------------------------

      subroutine day2nsec(cal_year,nsec,cal_month,cal_day,cal_hour,
     $     cal_min,cal_sec)
      use param
      use mod_mpi

!  calender day to nsec

      implicit none
      integer,intent(in) :: cal_year,cal_month,cal_day,
     &   cal_hour,cal_min,cal_sec
      integer,intent(out) :: nsec
      integer :: mday(12),msum,m

      mday(1)=31
      if(mod(cal_year,4)==0) then
         mday(2)=29
      else
         mday(2)=28
      endif
      mday(3)=31
      mday(4)=30
      mday(5)=31
      mday(6)=30
      mday(7:8)=31
      mday(9)=30
      mday(10)=31
      mday(11)=30
      mday(12)=31

      msum=0
      if(cal_month /= 1) then
      do m=1,cal_month-1
         msum=msum+mday(m)
      enddo
      endif

      msum = msum + cal_day-1

      nsec = msum*86400+cal_hour*3600+cal_min*60+cal_sec

      return
      end subroutine day2nsec

!------------------------------------------------------
      subroutine new_year
      use param
      use mod_mpi
#ifdef NCIO
      use netcdf
#endif
      implicit none
#include "common.h"

#ifdef NCEPSF
#ifndef EVAPO
      real(4) :: evapo_d(im,jm,0:nsf+1)
#endif
#endif

#ifdef NCIO
      integer :: ncstat
#endif
      integer :: nday_prev,i,j,n

      nsec=nsec-nday*86400

#ifdef TESTNY
      if(sfc_clim==1) return
      if(ip.eq.imaster) write(*,'("New Year: ",i4,"-->",i4)')
     $     cal_year,cal_year+1
      cal_year=cal_year+1
      nday_prev=nday
      if(mod(cal_year,4)==0) then
         nday=366
      else
         nday=365
      endif
      nxyear=nday*60*60*24/idnint(dtuv)
      nxmonth=nxyear/12
#ifdef NESTED
      if(obc_clim==1) then
         do n=1,nrbd
            ahourbd(n)=ahourbd(n)*dble(nday)/dble(nday_prev)
         enddo
      endif
#endif
!     ---   read next year's sfc
      do j=1,jml
      do i=1,iml
         wsx_d(i,j,0)=wsx_d(i,j,nday_prev)
         wsy_d(i,j,0)=wsy_d(i,j,nday_prev)
         sc_wind_d(i,j,0)=sc_wind_d(i,j,nday_prev)
         stdev_w_d(i,j,0)=stdev_w_d(i,j,nday_prev)
         dpt_t2m_d(i,j,0)=dpt_t2m_d(i,j,nday_prev)
         temp_2m_d(i,j,0)=temp_2m_d(i,j,nday_prev)
#ifdef NCEPSF
         dswrf_d(i,j,0)=dswrf_d(i,j,nday_prev)
#ifdef DLWRF
         dlwrf_d(i,j,0)=dlwrf_d(i,j,nday_prev)
#else
         cloud_d(i,j,0)=cloud_d(i,j,nday_prev)
#endif !end DLWRF
#else
         cloud_d(i,j,0)=cloud_d(i,j,nday_prev)
         tot_sol_d(i,j,0)=tot_sol_d(i,j,nday_prev)
#endif !end NCEPSF
         wflux_d(i,j,0)=wflux_d(i,j,nday_prev)
      enddo
      enddo
#ifdef NCIO
      if(ip.eq.imaster) then
         ncstat=nf90_open(sffile2,NF90_NOWRITE,bnddt)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
      endif
      call read_nc_fd(wsx_d(1,1,1),bnddt,1,nday,'wsx')
      call read_nc_fd(wsy_d(1,1,1),bnddt,1,nday,'wsy')
      call read_nc_fd(sc_wind_d(1,1,1),bnddt,1,nday,'sc_wind')
      call read_nc_fd(dpt_t2m_d(1,1,1),bnddt,1,nday,'dpt_t2m')
      call read_nc_fd(temp_2m_d(1,1,1),bnddt,1,nday,'temp_2m')
#ifdef NCEPSF
      call read_nc_fd(dswrf_d(1,1,1),bnddt,1,nday,'dswrf')
#ifdef DLWRF
      call read_nc_fd(dlwrf_d(1,1,1),bnddt,1,nday,'dlwrf')
#else
      call read_nc_fd(cloud_d(1,1,1),bnddt,1,nday,'cloud')
#endif !end DLWRF
#endif !end NCEPSF
      call read_nc_fd(wflux_d(1,1,1),bnddt,1,nday,'wflux')
#ifndef EVAPO
      call read_nc_fd(evapo_d(1,1,1),bnddt,1,nday,'evapo')
      do n=1,nday
      do j=1,jml
      do i=1,iml
         wflux_d(i,j,n)=(wflux_d(i,j,n)-evapo_d(i,j,n))*tex(i,j,1)
      enddo
      enddo
      enddo
#endif
      if(ip.eq.imaster) ncstat=nf90_close(bnddt)
#else ! not NCIO
      if(ip.eq.imaster) rewind(bnddt)
      call read_fd(wsx_d,bnddt)
      call read_fd(wsy_d,bnddt)
      call read_fd(sc_wind_d,bnddt)
      call read_fd(stdev_w_d,bnddt)
      call read_fd(dpt_t2m_d,bnddt)
      call read_fd(temp_2m_d,bnddt)
#ifdef NCEPSF
      call read_fd(dswrf_d,bnddt)
      call read_fd(dlwrf_d,bnddt)
#else
      call read_fd(cloud_d,bnddt)
      call read_fd(tot_sol_d,bnddt)
#endif !end NCEPSF
      call read_fd(wflux_d,bnddt)
      bnddt=bnddt+1
#endif
#endif !end TESTNY

      return
      end subroutine new_year

!--------------------------------------------------
      subroutine fac_monthly(cal_year,nsec,fca,fcb,mnt,mntp1)

! nsec to calender day

      use param
      use mod_mpi

      implicit none
      integer,intent(in) :: nsec,cal_year
      integer,intent(out) :: mnt,mntp1
      real*8,intent(out) :: fca,fcb
      integer :: mday(12),msum,m,cal_day,cal_month
      real*8 c_day


      mday(1)=31
      if(mod(cal_year,4)==0) then
         mday(2)=29
      else
         mday(2)=28
      endif
      mday(3)=31
      mday(4)=30
      mday(5)=31
      mday(6)=30
      mday(7:8)=31
      mday(9)=30
      mday(10)=31
      mday(11)=30
      mday(12)=31

      c_day=dble(nsec)/86400.
      cal_day=int(c_day)+1
      msum=1
      do m=1,12
         if(cal_day>=msum .and. cal_day<msum+mday(m))then
            cal_month=m
            exit
         endif
         msum=msum+mday(m)
      enddo
      c_day=c_day-dble(msum)+1.

      if(c_day*2 > mday(cal_month) ) then
        mnt = cal_month
        mntp1 = cal_month+1
        if(mntp1 == 13) mntp1 = 1
        fca = (2.*c_day-dble( mday(mnt) ))
     &    /dble( mday(mnt)+mday(mntp1) )
        fcb= 1.d0-fca
      else
        mnt = cal_month-1
        mntp1 = cal_month
        if(mnt == 0 ) mnt = 12
        fca = (2.*c_day+dble( mday(mnt) ))
     &    /dble( mday(mnt)+mday(mntp1) )
        fcb = 1.d0 -fca
      endif

      return
      end subroutine fac_monthly
