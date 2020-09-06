
! subroutine to make reference temperature and salinity
!    be interpolation from monthly data

!  $Id: ref_month.F90 5 2008-06-11 13:23:18Z ishikawa $

      subroutine ref_month

      use param
      use mod_mpi
      implicit none

#include "common.h"

      integer mnt,isday,mntp1
      integer i,j,k
      real*8 fca,fcb,aday

#ifdef GL11M
!  monthly forcing

      nnday = mod(nkai,360*3600*24/nint(dtuv))
      mnt = ((nnday-1)*nint(dtuv)/3600/24+15)/30
      isday = (mnt*30-15)*3600*24/nint(dtuv)
      fca = dble(nnday-isday)*dtuv/dble(3600*24*30)
      fcb = 1.-fca

      do k = 1,km
      do j = 1,jml
      do i = 1,iml
         tref(i,j,k)=tex(i,j,k)*(
     &        fcb*tref12(i,j,k,mnt)+fca*tref12(i,j,k,mnt+1))
         sref(i,j,k)=tex(i,j,k)*(
     &        fcb*sref12(i,j,k,mnt)+fca*sref12(i,j,k,mnt+1))
      enddo
      enddo
      enddo
#endif

#ifdef PC68M
!   monthly forcing

      if(nday==360) then
        nnday = mod(nkai,360*3600*24/nint(dtuv))
        mnt = ((nnday-1)*nint(dtuv)/3600/24+15)/30
        mntp1=mnt+1
        isday = (mnt*30-15)*3600*24/nint(dtuv)
        fca = dble(nnday-isday)*dtuv/dble(3600*24*30)
        fcb = 1.-fca
        if(mnt.eq.0)mnt=12
        if(mntp1.eq.13)mntp1=1
      else
        call fac_monthly(cal_year,nsec,fca,fcb,mnt,mntp1)
      endif

      do j = 1,jml
      do k = 1,km
      do i = 1,iml
         tref(i,j,k)=tex(i,j,k)*dble(
     &        fcb*tref12(i,j,k,mnt)+fca*tref12(i,j,k,mntp1))
         sref(i,j,k)=tex(i,j,k)*dble(
     &        fcb*sref12(i,j,k,mnt)+fca*sref12(i,j,k,mntp1))
      enddo
      enddo
      enddo
#endif

#ifdef ROKKA
      do j=1,jml
      do k=1,km
      do i=1,iml
         tref(i,j,k)=tex(i,j,k)*(tref12(i,j,k,mb1)
     $        +(tref12(i,j,k,mb2)-tref12(i,j,k,mb1))*cbf)
         sref(i,j,k)=tex(i,j,k)*(sref12(i,j,k,mb1)
     $        +(sref12(i,j,k,mb2)-sref12(i,j,k,mb1))*cbf)
      enddo
      enddo
      enddo
#endif
#ifdef NWNPAC
      do j=1,jml
      do k=1,km
      do i=1,iml
         tref(i,j,k)=tex(i,j,k)*(tref12(i,j,k,mb1)
     $        +(tref12(i,j,k,mb2)-tref12(i,j,k,mb1))*cbf)
         sref(i,j,k)=tex(i,j,k)*(sref12(i,j,k,mb1)
     $        +(sref12(i,j,k,mb2)-sref12(i,j,k,mb1))*cbf)
      enddo
      enddo
      enddo
#endif
#ifdef JP68M
      if(nday==360) then
        aday = dble(nsec)/3600./24.
        mnt = floor( (aday-15) / 30 ) + 1
        mntp1=mnt+1
        fca = (nsec - dble(mntp1*30-15)*3600.*24.)/(30.*24.*3600.)
        fcb = 1.-fca
        if(mnt == 0 ) mnt  = 12
        if(mntp1==13) mntp1= 1
      else
        call fac_monthly(cal_year,nsec,fca,fcb,mnt,mntp1)
      endif

      do j=1,jml
      do k=1,km
      do i=1,iml
         tref(i,j,k)=tex(i,j,k)*dble(
     &        fcb*tref12(i,j,k,mnt)+fca*tref12(i,j,k,mntp1))
         sref(i,j,k)=tex(i,j,k)*dble(
     &        fcb*sref12(i,j,k,mnt)+fca*sref12(i,j,k,mntp1))
      enddo
      enddo
      enddo
#endif

      return
      end
