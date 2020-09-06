ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine write                                              c
c
c  $Id: write.F90 11 2008-12-11 05:08:04Z ishikawa $
c                                                                   c
c     write real*4 data                                             c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write


      use param
      use mod_mpi
#ifdef NCIO
      use netcdf
#endif
      implicit none

#include "common.h"


      real*8 xw(im,jm,0:km+1),xhice(im,jm)
      real*8 xdepml(im,jm),sst,drdn,drup,dpth,aday,ayear
      integer i,j,k,iyear,imin,ihour,iday,month
#ifdef NCIO
      integer :: ncvar,ncstat
#endif
      do k=0,km+1
      do j=1,jml
      do i=1,iml
         xw(i,j,k)=wl(i,j,k)/(areat(i,j,k)+1.-tex(i,j,k))
      enddo
      enddo
      enddo

      do j = 1,jm
      do i = 1,im
        xdepml(i,j) = 0.
      enddo
      enddo

      do k = km,kmlmin,-1
      do j = 3,jml-2
      do i = 3,iml-2
         sst = dmax1(t(i,j,1),-3.d0)
         sst = dmin1(sst,30.d0)
         rhomix(i,j) = 1.d0-tex(i,j,1)+tex(i,j,1)*(
     &        rhoo(i,j,kmlmin)+drfit(1) +drfit(2)*sst
     &        +drfit(3)*sst*sst+drfit(4)*sst**3
     &        +drfit(5)*sst**4 )

         drdn = rhomix(i,j) - rhoo(i,j,k)
         drup = rhomix(i,j) - rhoo(i,j,k-1)
         if(drdn.lt.0.d0.and.drup.gt.0.d0)then
            if(k-1.le.kadp)then
               dpth=dep(k-1)*(1.d0+hcl(i,j)/depadp)+0.5d0*dzt(i,j,k-1)
            else
               dpth=dep(k-1)+hcl(i,j)+0.5d0*dzt(i,j,k-1)
            endif
            xdepml(i,j) = dpth + 0.5d0*(dzt(i,j,k-1)+dzt(i,j,k))
     &           *drup/(rhoo(i,j,k)-rhoo(i,j,k-1))
        endif

        xdepml(i,j)=dmax1(xdepml(i,j),10.d2)
        xdepml(i,j)=dmin1(xdepml(i,j),2000.d2)*tex(i,j,1)
      enddo
      enddo
      enddo

      do j=1,jml
      do i=1,iml
         if(aice(i,j)*areat(i,j,1).gt.0.)then
            xhice(i,j)=volice(i,j)/(aice(i,j)*areat(i,j,1))
         else
            xhice(i,j)=0.
         endif
      enddo
      enddo
c
c
#ifdef PC68M
	  rewind(continu)
      if(ip.eq.imaster) then
        write(result) ahour
      endif
      call write_3d8_2d4_r4(u,result)
      call write_3d8_2d4_r4(v,result)
      call write_3d8_2d4_r4(xw,result)
#ifndef T3OUT
      call write_3d8_2d4_r4(t,result)
      call write_3d8_2d4_r4(s,result)
      call write_3d8_2d4_r4(rhoo,result)
      call write_2d8_r4(hcl,result)
      call write_2d8_r4(sfu,result)
      call write_2d8_r4(sfv,result)
      call write_3d8_2d4_r4(tke,result)
c      call write_3d8_2d4_r4(vddt,result)
c      call write_3d8_2d4_r4(vdds,result)
      call write_2d8_r4(xdepml,result)
      call write_2d8_r4(aice,result)
      call write_2d8_r4(xhice,result)
      call write_2d8_r4(uice,result)
      call write_2d8_r4(vice,result)
#endif
c
      result=result+1
#else !PC68M
#ifdef NCIO
      ncnt_rs=ncnt_rs+1

      if(ip==imaster) then

         ncstat=nf90_inq_varid(result,'ahour',ncvar)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'ahour'
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

         ncstat=nf90_put_var(result,ncvar,ahour,ncnt_rs)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

#ifndef CLIMAT
        ncstat=nf90_inq_varid(result,'nsec',ncvar)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nsec'
         write(*,*) nf90_strerror(ncstat)
      endif
#endif
         ncstat=nf90_put_var(result,ncvar,nsec,ncnt_rs)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'nsec'
            write(*,*) nf90_strerror(ncstat)
         endif
#endif
#endif !CLIMAT

      endif
      call write_nc_3d8_r4(u,result,ncnt_rs,'u')
      call write_nc_3d8_r4(v,result,ncnt_rs,'v')
      call write_nc_3d8_r4(xw,result,ncnt_rs,'w')
#ifndef T3OUT
      call write_nc_3d8_r4(t,result,ncnt_rs,'t')
      call write_nc_3d8_r4(s,result,ncnt_rs,'s')
      call write_nc_3d8_r4(rhoo,result,ncnt_rs,'rhoo')
      call write_nc_2d8_r4(hcl,result,ncnt_rs,'hcl')
      call write_nc_2d8_r4(sfu,result,ncnt_rs,'sfu')
      call write_nc_2d8_r4(sfv,result,ncnt_rs,'sfv')
      call write_nc_3d8_r4(tke,result,ncnt_rs,'tke')
      call write_nc_2d8_r4(xdepml,result,ncnt_rs,'depml')
#ifdef ICE
      call write_nc_2d8_r4(aice,result,ncnt_rs,'aice')
      call write_nc_2d8_r4(xhice,result,ncnt_rs,'hice')
      call write_nc_2d8_r4(uice,result,ncnt_rs,'uice')
      call write_nc_2d8_r4(vice,result,ncnt_rs,'vice')
#endif
#endif

      if(ip.eq.imaster) then
         ncstat=nf90_sync(result)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif
      endif
#else !NCIO
      if(ip.eq.imaster) then
        write(result) ahour
      endif
      call write_3d8_r4(u,result)
      call write_3d8_r4(v,result)
      call write_3d8_r4(xw,result)
#ifndef T3OUT
      call write_3d8_r4(t,result)
      call write_3d8_r4(s,result)
      call write_3d8_r4(rhoo,result)
      call write_2d8_r4(hcl,result)
      call write_2d8_r4(sfu,result)
      call write_2d8_r4(sfv,result)
      call write_3d8_r4(tke,result)
c      call write_3d8_r4(vddt,result)
c      call write_3d8_r4(vdds,result)
      call write_2d8_r4(xdepml,result)
#ifdef ICE
      call write_2d8_r4(aice,result)
      call write_2d8_r4(xhice,result)
      call write_2d8_r4(uice,result)
      call write_2d8_r4(vice,result)
#endif
#endif
#endif !NCIO
#endif !PC68
c
c
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
!      imin=(ahour-iyear*8640.d0)*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
        write(6,*)nkai,'step : written on result file'
        write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &       ' hour=',ihour,' min=',imin
        write(6,*)' '
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine wrcontin                                           c
c
c     save for continue                                             c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wrcontin

      use param
      use mod_mpi

      implicit none
#include "common.h"

      real*8 aday,ayear
      integer iyear,month,iday,ihour,imin,isec

      nkai=nkai-1
c
#ifdef PC68M
      if(ip.eq.imaster) then
         rewind(continu)
!         write(continu) nkai,ahour
      if(nsec>=nday*86400) nsec=nsec-nday*86400
         write(continu) nkai,ahour,nsec
#ifdef DEBUG
        write(*,*) 'last step: ',nkai,ahour,nsec
#endif
         write(continu) pd,pm,ddmna,dmn
      endif
      call write_3d8_2d8_r8(u,continu)
      call write_3d8_2d8_r8(v,continu)
      call write_3d8_2d8_r8(t,continu)
      call write_3d8_2d8_r8(s,continu)
      call write_3d8_2d8_r8(tke,continu)
      call write_2d8_r8(hcl,continu)
      call write_2d8_r8(um,continu)
      call write_2d8_r8(vm,continu)
      call write_2d8_r8(aice,continu)
      call write_2d8_r8(volice,continu)
      call write_2d8_r8(uice,continu)
      call write_2d8_r8(vice,continu)
#else
      if(ip.eq.imaster) then
#ifdef NESTED
      rewind(continu)
!        ctime=ctime-dtts
      if(nsec>=nday*86400) nsec=nsec-nday*86400
      write(continu) nkai,ahour,nsec
!         write(continu) nkai,ahour,ctime_ini
#else
         write(continu) nkai,ahour
#endif
         write(continu) pd,pm,ddmna,dmn
      endif
      call write_3d8_r8(u,continu)
      call write_3d8_r8(v,continu)
      call write_3d8_r8(t,continu)
      call write_3d8_r8(s,continu)
      call write_3d8_r8(tke,continu)
      call write_2d8_r8(hcl,continu)
      call write_2d8_r8(um,continu)
      call write_2d8_r8(vm,continu)
      call write_2d8_r8(aice,continu)
      call write_2d8_r8(volice,continu)
      call write_2d8_r8(uice,continu)
      call write_2d8_r8(vice,continu)
#endif
c
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
         write(*,*)nkai,' step : written on continue file'
         write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &        ' hour=',ihour,' min=',imin
         write(*,*)' '
      endif

#ifdef NESTED
#ifndef CLIMAT
      call calender(cal_year,nsec,month,iday,ihour,imin,isec)
      if(ip==imaster) then
         write(*,'("date ",i4,i2.2,i2.2)') cal_year,month,iday
         write(*,'("time ",i2.2,":",i2.2," ",i2.2)') ihour,imin,isec
         write(*,*)
      endif
#endif
#endif

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                   c
!     subroutine wrmean                                             c
!                                                                   c
!     write mean data                                               c
!                                                                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wrmean

      use param
      use mod_mpi
#ifdef NCIO
      use netcdf
#endif

      implicit none
#include "common.h"

      real*8 xdepml(im,jm),sst,drdn,drup,dpth,aday,ayear
      integer iyear,month,iday,ihour,imin,i,j,k
#ifdef NCIO
      integer :: ncstat,ncvar
#endif
#ifndef CLIMAT
      integer,save :: nkai_rm=0
#endif
#if defined(ROKKA) || defined(NWNPAC)
      real(8),save :: wmnsec=0.d0
#endif
#if defined(ROKKA) && defined(EXP2007)
      integer,save :: ncnt_fm=0,nfile=1
      integer,parameter :: ndpf=8
#endif
#ifdef CLIAMT
      if(mod(nkai,writ_mn).eq.1)then
#else
      nkai_rm=nkai_rm+1
      if(mod(nkai_rm,writ_mn)==1) then
#endif

         wmahour=0.
         wmnsec=0.d0

         do k=0,km+1
         do j=1,jml
         do i=1,iml
            wmu(i,j,k)=0.
            wmv(i,j,k)=0.
            wmt(i,j,k)=0.
            wms(i,j,k)=0.
            wmrho(i,j,k)=0.
            wmwl(i,j,k)=0.
            wmtke(i,j,k)=0.
            wmvdts(i,j,k)=0.
            wmvddt(i,j,k)=0.
            wmvdds(i,j,k)=0.
         enddo
         enddo
         enddo

         do j=1,jml
         do i=1,iml
            wmhcl(i,j)=0.
            wmsfu(i,j)=0.
            wmsfv(i,j)=0.
            wmmld(i,j)=0.
            wmaice(i,j)=0.
            wmhice(i,j)=0.
            wmuice(i,j)=0.
            wmvice(i,j)=0.
            wmheat(i,j)=0.
            wmwater(i,j)=0.
            wmheat_r(i,j)=0.
            wmsal_r(i,j)=0.
         enddo
         enddo

      endif

      wmahour=wmahour+ahour/dble(writ_mn)
#if defined(ROKKA) || defined(NWNPAC)
      wmnsec=wmnsec+dble(nsec)/dble(writ_mn)
#endif
c
      do j = 3,jml-2
      do k = km,kmlmin,-1
      do i = 3,iml-2
         sst = dmax1(t(i,j,1),-3.d0)
         sst = dmin1(sst,30.d0)
         rhomix(i,j) = 1.d0-tex(i,j,1)+tex(i,j,1)*(
     &        rhoo(i,j,kmlmin)+drfit(1) +drfit(2)*sst
     &        +drfit(3)*sst*sst+drfit(4)*sst**3
     &        +drfit(5)*sst**4 )

         drdn = rhomix(i,j) - rhoo(i,j,k)
         drup = rhomix(i,j) - rhoo(i,j,k-1)
         if(drdn.lt.0.d0.and.drup.gt.0.d0)then
            if(k-1.le.kadp)then
               dpth=dep(k-1)*(1.d0+hcl(i,j)/depadp)+0.5d0*dzt(i,j,k-1)
            else
               dpth=dep(k-1)+hcl(i,j)+0.5d0*dzt(i,j,k-1)
            endif
            xdepml(i,j) = dpth + 0.5d0*(dzt(i,j,k-1)+dzt(i,j,k))
     &           *drup/(rhoo(i,j,k)-rhoo(i,j,k-1))
        endif

        xdepml(i,j)=dmax1(xdepml(i,j),10.d2)
        xdepml(i,j)=dmin1(xdepml(i,j),2000.d2)*tex(i,j,1)
      enddo
      enddo
      enddo

      do k=0,km+1
      do j=1,jml
      do i=1,iml
         wmu(i,j,k)=wmu(i,j,k)+ub(i,j,k)/real(writ_mn)
         wmv(i,j,k)=wmv(i,j,k)+vb(i,j,k)/real(writ_mn)
         wmt(i,j,k)=wmt(i,j,k)+tb(i,j,k)/real(writ_mn)
         wms(i,j,k)=wms(i,j,k)+sb(i,j,k)/real(writ_mn)
         wmrho(i,j,k)=wmrho(i,j,k)+rhoo(i,j,k)/real(writ_mn)
         wmwl(i,j,k)=wmwl(i,j,k)
     &        +wl(i,j,k)/(areat(i,j,k)+1.-tex(i,j,k))/real(writ_mn)
         wmtke(i,j,k)=wmtke(i,j,k)+tkeb(i,j,k)/real(writ_mn)
         wmvdts(i,j,k)=wmvdts(i,j,k)+vdts(i,j,k)/real(writ_mn)
         wmvddt(i,j,k)=wmvddt(i,j,k)+vddt(i,j,k)/real(writ_mn)
         wmvdds(i,j,k)=wmvdds(i,j,k)+vdds(i,j,k)/real(writ_mn)
      enddo
      enddo
      enddo

      do j=1,jml
      do i=1,iml
         wmhcl(i,j)=wmhcl(i,j)+hcl(i,j)/real(writ_mn)
         wmsfu(i,j)=wmsfu(i,j)+sfund(i,j)/real(writ_mn)
         wmsfv(i,j)=wmsfv(i,j)+sfvnd(i,j)/real(writ_mn)
         wmmld(i,j)=wmmld(i,j)+xdepml(i,j)/real(writ_mn)
         wmaice(i,j)=wmaice(i,j)+aiceb(i,j)/real(writ_mn)
         if(aiceb(i,j).gt.0.d0)then
             wmhice(i,j)=wmhice(i,j)+voliceb(i,j)
     &    /(areat(i,j,1)+tex(i,j,1)-1.)/aiceb(i,j)/real(writ_mn)
         endif
         wmuice(i,j)=wmuice(i,j)+uice(i,j)/real(writ_mn)
         wmvice(i,j)=wmvice(i,j)+vice(i,j)/real(writ_mn)
         wmheat(i,j)=wmheat(i,j)
     &        +tex(i,j,1)*(netq_o(i,j)+netq_i(i,j)
!     &        +gref(i,j,1)*(tref(i,j,1)-t(i,j,1))
!     &        *dzt(i,j,1)*ddmna/0.24d0
     &        )/real(writ_mn)
         wmwater(i,j)=wmwater(i,j)+wflux(i,j)/real(writ_mn)
         wmheat_r(i,j)=wmheat_r(i,j)
     &        +tex(i,j,1)*gref(i,j,1)*(tref(i,j,1)-t(i,j,1))
     &        *dzt(i,j,1)*ddmna/0.24d0/real(writ_mn)
         wmsal_r(i,j)=wmsal_r(i,j)
     &        +tex(i,j,1)*gref(i,j,1)*(sref(i,j,1)-s(i,j,1))
     &        *dzt(i,j,1)/real(writ_mn)
      enddo
      enddo

#ifdef CLIMAT
      if(mod(nkai,writ_mn).eq.0)then
#else
      if(mod(nkai_rm,writ_mn).eq.0)then
#endif

#ifdef NCIO
      ncnt_rm=ncnt_rm+1

#if defined(ROKKA) && defined(EXP2007)


      ncnt_fm=ncnt_fm+1

      if(ip==imaster) then

         ncstat=nf90_inq_varid(resmean(nfile),'nsec',ncvar)

         ncstat=nf90_put_var(resmean(nfile),ncvar,nint(wmnsec),ncnt_rm)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'nsec',resmean(nfile),nfile
            write(*,*) nf90_strerror(ncstat)
!            stop
         endif
#endif
      endif
      call write_nc_3d4_r4(wmu,resmean(nfile),ncnt_rm,'u')
      call write_nc_3d4_r4(wmv,resmean(nfile),ncnt_rm,'v')
      call write_nc_3d4_r4(wmwl,resmean(nfile),ncnt_rm,'w')
#ifndef BNDTIDE_MAINDATA
      call write_nc_3d4_r4(wmt,resmean(nfile),ncnt_rm,'t')
      call write_nc_3d4_r4(wms,resmean(nfile),ncnt_rm,'s')
#endif
      call write_nc_3d4_r4(wmrho,resmean(nfile),ncnt_rm,'rhoo')
      call write_nc_2d4_r4(wmhcl,resmean(nfile),ncnt_rm,'hcl')
      call write_nc_2d4_r4(wmsfu,resmean(nfile),ncnt_rm,'sfu')
      call write_nc_2d4_r4(wmsfv,resmean(nfile),ncnt_rm,'sfv')
#ifndef BNDTIDE_MAINDATA
      call write_nc_3d4_r4(wmtke,resmean(nfile),ncnt_rm,'tke')
      call write_nc_2d4_r4(wmmld,resmean(nfile),ncnt_rm,'depml')
      call write_nc_3d4_r4(wmvdts,resmean(nfile),ncnt_rm,'vdts')
#endif

      if(ip==imaster) then
         ncstat=nf90_sync(resmean(nfile))
#ifdef DEBUG
         if(ncstat/=0) write(*,*) nf90_strerror(ncstat)
#endif
      endif

      if(ncnt_rm(1) >= ndpf) then
         if(ip==imaster) then
            ncstat=nf90_close(resmean(nfile))
         endif
         ncnt_rm=0
         nfile=nfile+1
      endif

      call write_nc_2d4_r4(wmheat,flxmean,ncnt_fm,'heat')
      call write_nc_2d4_r4(wmwater,flxmean,ncnt_fm,'fwater')


#else ! defined(ROKKA) && defined(EXP2007)
      if(ip.eq.imaster) then
         ncstat=nf90_sync(resmean)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

         ncstat=nf90_inq_varid(resmean,'ahour',ncvar)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'ahour'
            write(*,*) nf90_strerror(ncstat)
!            stop
         endif
#endif

         ncstat=nf90_put_var(resmean,ncvar,wmahour,ncnt_rm)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
!            stop
         endif
#endif

#if defined(ROKKA) || defined(NWNPAC)

         ncstat=nf90_inq_varid(resmean,'nsec',ncvar)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'nsec'
            write(*,*) nf90_strerror(ncstat)
!            stop
         endif
#endif

         ncstat=nf90_put_var(resmean,ncvar,nint(wmnsec),ncnt_rm)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'nsec'
            write(*,*) nf90_strerror(ncstat)
!            stop
         endif
#endif

#endif !end defined(ROKKA) || defined(NWNPAC)
      endif
      call write_nc_3d4_r4(wmu,resmean,ncnt_rm,'u')
      call write_nc_3d4_r4(wmv,resmean,ncnt_rm,'v')
      call write_nc_3d4_r4(wmwl,resmean,ncnt_rm,'w')
#ifndef BNDTIDE_MAINDATA
      call write_nc_3d4_r4(wmt,resmean,ncnt_rm,'t')
      call write_nc_3d4_r4(wms,resmean,ncnt_rm,'s')
#endif
      call write_nc_3d4_r4(wmrho,resmean,ncnt_rm,'rhoo')
      call write_nc_2d4_r4(wmhcl,resmean,ncnt_rm,'hcl')
      call write_nc_2d4_r4(wmsfu,resmean,ncnt_rm,'sfu')
      call write_nc_2d4_r4(wmsfv,resmean,ncnt_rm,'sfv')
#ifndef BNDTIDE_MAINDATA
      call write_nc_3d4_r4(wmtke,resmean,ncnt_rm,'tke')
      call write_nc_2d4_r4(wmmld,resmean,ncnt_rm,'depml')
#if defined(SGOUT) && defined(EXP2006)
      call write_nc_3d4_r4(wmvdts,resmean,ncnt_rm,'vdts')
#endif
#ifdef ICE
      call write_nc_2d4_r4(wmaice,resmean,ncnt_rm,'aice')
      call write_nc_2d4_r4(wmhice,resmean,ncnt_rm,'hice')
      call write_nc_2d4_r4(wmuice,resmean,ncnt_rm,'uice')
      call write_nc_2d4_r4(wmvice,resmean,ncnt_rm,'vice')
#endif
#endif
      call write_nc_2d4_r4(wmheat,flxmean,ncnt_rm,'heat')
      call write_nc_2d4_r4(wmwater,flxmean,ncnt_rm,'fwater')
#ifdef DEBUG11
      if(ip.eq.imaster) ncstat=nf90_sync(resmean)
#endif
#endif ! end not defined(ROKKA) && defined(EXP2007)

#else ! not NCIO
#ifdef PC68M
      if(ip.eq.imaster) then
        write(resmean)wmahour
      endif
      call write_3d4_2d4_r4(wmu,resmean)
      call write_3d4_2d4_r4(wmv,resmean)
      call write_3d4_2d4_r4(wmwl,resmean)
      call write_3d4_2d4_r4(wmt,resmean)
      call write_3d4_2d4_r4(wms,resmean)
      call write_3d4_2d4_r4(wmrho,resmean)
      call write_2d4_r4(wmhcl,resmean)
      call write_2d4_r4(wmsfu,resmean)
      call write_2d4_r4(wmsfv,resmean)
      call write_3d4_2d4_r4(wmtke,resmean)
!      call write_3d4_2d4_r4(wmvddt,resmean)
!      call write_3d4_2d4_r4(wmvdds,resmean)
      call write_2d4_r4(wmmld,resmean)
      call write_2d4_r4(wmaice,resmean)
      call write_2d4_r4(wmhice,resmean)
      call write_2d4_r4(wmuice,resmean)
      call write_2d4_r4(wmvice,resmean)

      resmean=resmean+1
#else !PC68M
#ifdef GL11M
      if(ip.eq.imaster) then
        write(resmean)wmahour
      endif
      call write_3d4_r4(wmu,resmean)
      call write_3d4_r4(wmv,resmean)
      call write_3d4_r4(wmwl,resmean)
      call write_3d4_r4(wmt,resmean)
      call write_3d4_r4(wms,resmean)
      call write_3d4_r4(wmrho,resmean)
      call write_2d4_r4(wmhcl,resmean)
      call write_2d4_r4(wmsfu,resmean)
      call write_2d4_r4(wmsfv,resmean)
      call write_3d4_r4(wmtke,resmean)
!      call write_3d4_r4(wmvddt,resmean)
!      call write_3d4_r4(wmvdds,resmean)
      call write_2d4_r4(wmmld,resmean)
      call write_2d4_r4(wmaice,resmean)
      call write_2d4_r4(wmhice,resmean)
      call write_2d4_r4(wmuice,resmean)
      call write_2d4_r4(wmvice,resmean)

      call write_2d4_r4(wmheat,flxmean)
      call write_2d4_r4(wmwater,flxmean)
      call write_2d4_r4(wmheat_r,flxmean)
      call write_2d4_r4(wmsal_r,flxmean)

#else ! GL11M
      if(ip.eq.imaster) then
        write(resmean)wmahour
      endif
      call write_3d4_r4(wmu,resmean)
      call write_3d4_r4(wmv,resmean)
      call write_3d4_r4(wmwl,resmean)
      call write_3d4_r4(wmt,resmean)
      call write_3d4_r4(wms,resmean)
      call write_3d4_r4(wmrho,resmean)
      call write_2d4_r4(wmhcl,resmean)
      call write_2d4_r4(wmsfu,resmean)
      call write_2d4_r4(wmsfv,resmean)
      call write_3d4_r4(wmtke,resmean)
      call write_3d4_r4(wmvddt,resmean)
      call write_3d4_r4(wmvdds,resmean)
      call write_2d4_r4(wmmld,resmean)
#ifdef ICE
      call write_2d4_r4(wmaice,resmean)
      call write_2d4_r4(wmhice,resmean)
      call write_2d4_r4(wmuice,resmean)
      call write_2d4_r4(wmvice,resmean)
#endif
#endif ! end not GL11M
#endif ! end not PC68M

#endif !end not NCIO

      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
!      imin=(ahour-iyear*8640.d0)*60.d0+.01d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
         write(*,*)nkai,' step : written on mean-state file'
         write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &        ' hour=',ihour,' min=',imin
         write(*,*)' '
      endif

      endif
      return
      end
#ifdef BDOUT
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine odwrite                                            c
c                                                                   c
c     write open boundary data                                      c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine odwrite

      use param
      use mod_mpi

      implicit none
#include "common.h"
      integer,save :: ncount_bd=0
      integer :: i,j,k,imin,ihour,iday,month,iyear
      real(8) :: aday,ayear
c
c
!      if(mod(nkai,nwrod).eq.1)then
      ncount_bd=ncount_bd+1
!      if(ip==imaster) write(*,*) 'ncount_bd=',ncount_bd
      if(mod(ncount_bd,nwrod)==1) then
c
      ahour_od=0.d0
      ddmnar_od=0.d0
c
      do k=0,km+1
         pd_od(k)=0.d0
      enddo
c
      do k=0,km+1
      do j=1,jml
      do i=1,iml
c         u_od(i,j,k)=0.d0
c         v_od(i,j,k)=0.d0
         t_od(i,j,k)=0.d0
         s_od(i,j,k)=0.d0
         tke_od(i,j,k)=0.d0
      enddo
      enddo
      enddo
c
      do j=1,jml
      do i=1,iml
         hcl_od(i,j)=0.d0
         um_od(i,j)=0.d0
         vm_od(i,j)=0.d0
         aice_od(i,j)=0.d0
         volice_od(i,j)=0.d0
         uice_od(i,j)=0.d0
         vice_od(i,j)=0.d0
      enddo
      enddo
c
      endif
c
      ahour_od=ahour_od+ahour/dble(nwrod)
      ddmnar_od=ddmnar_od+ddmnar/dble(nwrod)
c
      do k=0,km+1
         pd_od(k)=pd_od(k)+pd(k)/dble(nwrod)
      enddo
c
      do k=0,km+1
      do j=1,jml
      do i=1,iml
c         u_od(i,j,k)=u_od(i,j,k)+ub(i,j,k)/dble(nwrod)
c         v_od(i,j,k)=v_od(i,j,k)+vb(i,j,k)/dble(nwrod)
         t_od(i,j,k)=t_od(i,j,k)+tb(i,j,k)/dble(nwrod)
         s_od(i,j,k)=s_od(i,j,k)+sb(i,j,k)/dble(nwrod)
         tke_od(i,j,k)=tke_od(i,j,k)+tkeb(i,j,k)/dble(nwrod)
      enddo
      enddo
      enddo
c
      do j=1,jml
      do i=1,iml
         hcl_od(i,j)=hcl_od(i,j)+hcl(i,j)/dble(nwrod)
         um_od(i,j)=um_od(i,j)+umb(i,j)/dble(nwrod)
         vm_od(i,j)=vm_od(i,j)+vmb(i,j)/dble(nwrod)
         aice_od(i,j)=aice_od(i,j)+aiceb(i,j)/dble(nwrod)
         volice_od(i,j)=volice_od(i,j)+voliceb(i,j)/dble(nwrod)
         uice_od(i,j)=uice_od(i,j)+uice(i,j)/dble(nwrod)
         vice_od(i,j)=vice_od(i,j)+vice(i,j)/dble(nwrod)
      enddo
      enddo
c
!      if(mod(nkai,nwrod).eq.0)then
      if(mod(ncount_bd,nwrod).eq.0)then
c
      if(ip.eq.imaster) then
         write(fuoutdt) ahour_od
         write(fuoutdt) pd_od,ddmnar_od

      endif
c      call write_3d8_r8(u_od,fuoutdt)
c      call write_3d8_r8(v_od,fuoutdt)
      call write_3d8_r8(t_od,fuoutdt)
      call write_3d8_r8(s_od,fuoutdt)
      call write_3d8_r8(tke_od,fuoutdt)
      call write_2d8_r8(hcl_od,fuoutdt)
      call write_2d8_r8(um_od,fuoutdt)
      call write_2d8_r8(vm_od,fuoutdt)
      call write_2d8_r8(aice_od,fuoutdt)
      call write_2d8_r8(volice_od,fuoutdt)
      call write_2d8_r8(uice_od,fuoutdt)
      call write_2d8_r8(vice_od,fuoutdt)

c
#if defined(NWNPAC) && defined(EXP2006)
      if(nsec>=nday*86400) nsec=nsec-nday*86400
      if(ip.eq.imaster) then
         write(outini) nkai,ahour,nsec
         write(outini) pd,pm,ddmna,dmn
      endif
      call write_3d8_r8(u,outini)
      call write_3d8_r8(v,outini)
      call write_3d8_r8(t,outini)
      call write_3d8_r8(s,outini)
      call write_3d8_r8(tke,outini)
      call write_2d8_r8(hcl,outini)
      call write_2d8_r8(um,outini)
      call write_2d8_r8(vm,outini)
      call write_2d8_r8(aice,outini)
      call write_2d8_r8(volice,outini)
      call write_2d8_r8(uice,outini)
      call write_2d8_r8(vice,outini)
#endif

c
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
!      imin=(ahour-iyear*8640.d0)*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30
c
      if(ip.eq.imaster) then
         write(*,*)nkai,' step : written on open bounday file'
         write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &        ' hour=',ihour,' min=',imin
         write(*,*)' '
      endif
c
      endif
c
      return
      end
c

#endif

#if defined(SGOUT) && !defined(EXP2006)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine sgwrite                                            c
c                                                                   c
c     write Sea-Gearn data                                          c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sgwrite

      use param
      use mod_mpi

      implicit none
#include "common.h"
      integer,save :: ncount_sg=0

      integer :: i,j,k,imin,ihour,iday,month,iyear
      real(8) :: aday,ayear

      ncount_sg=ncount_sg+1

      if(mod(ncount_sg,nwrsg).eq.1)then
!      if(mod(nkai,nwrsg).eq.1)then
c
      ahour_sg=0.d0
c
      do k=0,km+1
      do j=1,jml
      do i=1,iml
         u_sg(i,j,k)=0.d0
         v_sg(i,j,k)=0.d0
         w_sg(i,j,k)=0.d0
         vdts_sg(i,j,k)=0.d0
      enddo
      enddo
      enddo
c
      endif
c
      ahour_sg=ahour_sg+ahour/dble(nwrsg)
c
      do k=0,km+1
      do j=1,jml
      do i=1,iml
         u_sg(i,j,k)=u_sg(i,j,k)+ub(i,j,k)/dble(nwrsg)
         v_sg(i,j,k)=v_sg(i,j,k)+vb(i,j,k)/dble(nwrsg)
         w_sg(i,j,k)=w_sg(i,j,k)+wl(i,j,k)/dble(nwrsg)
         vdts_sg(i,j,k)=vdts_sg(i,j,k)+vdts(i,j,k)/dble(nwrsg)
      enddo
      enddo
      enddo
c
      if(mod(ncount_sg,nwrsg).eq.0)then
!      if(mod(nkai,nwrsg).eq.0)then
c
      do k=0,km+1
      do j=1,jml
      do i=1,iml
         w_sg(i,j,k)=w_sg(i,j,k)/(areat(i,j,k)+1.-tex(i,j,k))
      enddo
      enddo
      enddo
c
      if(ip.eq.imaster)then
         write(fuoutsg) ahour_sg
      endif
      call write_3d8_r4(u_sg,fuoutsg)
      call write_3d8_r4(v_sg,fuoutsg)
      call write_3d8_r4(w_sg,fuoutsg)
      call write_3d8_r4(vdts_sg,fuoutsg)
c

      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
!      imin=(ahour-iyear*8640.d0)*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30
c
      if(ip.eq.imaster) then
         write(*,*)nkai,' step : written on Sea-Gearn file'
         write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &        ' hour=',ihour,' min=',imin
         write(*,*)' '
      endif
c
      endif
c
      return
      end
c

#endif
!----------------------------------------------------
      subroutine wr_debug


      use param
      use mod_mpi

      implicit none
#include "common.h"

      real x3d(20,20,20),x2d(20,20)
      integer i,j,k,ix,jy

      jy = 115
      if(ip.eq.18) then
        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = u(i,j+jy,k)
        enddo
        enddo
        enddo
        write(75) x3d

        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = v(i,j+jy,k)
        enddo
        enddo
        enddo
        write(75) x3d

        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = t(i,j+jy,k)
        enddo
        enddo
        enddo
        write(75) x3d

        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = s(i,j+jy,k)
        enddo
        enddo
        enddo
        write(75) x3d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = hcl(i,j+jy)
        enddo
        enddo
        write(75) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = aice(i,j+jy)
        enddo
        enddo
        write(75) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = volice(i,j+jy)
        enddo
        enddo
        write(75) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = uice(i,j+jy)
        enddo
        enddo
        write(75) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = vice(i,j+jy)
        enddo
        enddo
        write(75) x2d
      endif

      ix = im-20
      if(ip.eq.17) then
        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = u(i+ix,j+jy,k)
        enddo
        enddo
        enddo
        write(76) x3d

        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = v(i+ix,j+jy,k)
        enddo
        enddo
        enddo
        write(76) x3d

        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = t(i+ix,j+jy,k)
        enddo
        enddo
        enddo
        write(76) x3d

        do k = 1,20
        do j = 1,20
        do i = 1,20
          x3d(i,j,k) = s(i+ix,j+jy,k)
        enddo
        enddo
        enddo
        write(76) x3d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = hcl(i+ix,j+jy)
        enddo
        enddo
        write(76) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = aice(i+ix,j+jy)
        enddo
        enddo
        write(76) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = volice(i+ix,j+jy)
        enddo
        enddo
        write(76) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = uice(i+ix,j+jy)
        enddo
        enddo
        write(76) x2d

        do j = 1,20
        do i = 1,20
          x2d(i,j) = vice(i+ix,j+jy)
        enddo
        enddo
        write(76) x2d

      endif

      return
      end
