ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine meanfc                                             c
c
c   $Id: mmchk.F90 2 2007-11-15 13:07:08Z ishikawa $
c                                                                   c
c     Diagnosis of model state.                                     c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mmchk

      use param
      use mod_mpi
      implicit none

#include "common.h"
c
      real*8 totarea,htmn,htmnpice,spm(km),eng,eng2
      integer iiday(km),i,j,k,iyear,month,iday,ihour,imin
      real*8 aday,ayear
      real*8 umax,vmax,tmax,smax,rmax,wmax,emax,
     &     umin,vmin,tmin,smin,wmin,emin,dtmin,dumin,
     &     aimax,volimax,himax,uimax,vimax,uimin,vimin,
     &     hclmax,hclmin,ummax,ummin,vmmax,vmmin,
     &     vdtsmax,vdtsmin,w
      integer iumax,jumax,iumin,jumin,ivmax,jvmax,ivmin,jvmin,
     &    itmax,jtmax,itmin,jtmin,ismax,jsmax,ismin,jsmin,
     &    iwmax,jwmax,iwmin,jwmin,iemax,jemax,iemin,jemin,
     &    idtmax,jdtmax,idtmin,jdtmin,iaimax,jaimax,iaimin,jaimin,
     &    ivolimax,jvolimax,ihimax,jhimax,
     &    iuimax,juimax,ivimax,jvimax,iuimin,juimin,ivimin,jvimin,
     &    ihclmax,jhclmax,ihclmin,jhclmin,iummax,jummax,iummin,jummin,
     &    ivmmax,jvmmax,ivmmin,jvmmin

      real*8 amax(km*7+8),amin(km*7+8)
      integer iamax(km*7+8),jamax(km*7+8),
     &     iamin(km*7+8),jamin(km*7+8)

!      real*8 test1(km),test2(km)

      if(ip.eq.imaster) then
        write(*,*) nkai,'step : max and min values'
      endif


      engina=0.d0
      engpta=0.d0
      enstro=0.d0
      ttmna=0.d0
      ssmna=0.d0
      do k=1,km
        tmn(k)=0.d0
        smn(k)=0.d0
        d2mn(k)=0.d0
        eni(k)=0.d0
        est(k)=0.d0
        htmn=0.d0
        spm(k)=0.d0
        dmn(k) = 0.d0
!        test1(k) = 0.
!        test2(k) = 0.
      enddo
c
      engexa=0.d0
      totarea=0.d0
      htmn = 0.d0
      htmnpice=0.d0
      do j=3,jml-2
      do i=3,iml-2
        engexa=engexa+.5d0*(sfun(i,j)*sfun(i,j)+sfvn(i,j)*sfvn(i,j))
     1    *areauu(j)*(hrr(i,j)+hclu(i,j))
        htmn=htmn+hcl(i,j)*areat(i,j,1)
        totarea=totarea+areat(i,j,1)
        htmnpice=htmnpice+hcl(i,j)*areat(i,j,1)
     &       +volice(i,j)*rho_i/rho_w
      enddo
      enddo
c
      call sum_mst1(engexa,1)
      call sum_mst1(htmn,1)
      call sum_mst1(totarea,1)
      call sum_mst1(htmnpice,1)
c
      if(ip.eq.imaster) then
        htmn=htmn/totarea
        htmnpice=htmnpice/totarea
        engexa=engexa/ttvol
      endif
c
      do j =2,jml
      do i =2,iml
        dzu(i,j,1)=ex(i,j,1)*(dz(1)+hclu(i,j))
        volt(i,j,1) = tex(i,j,1)*dz(1)
     $    *((tex(i-1,j-1,1)+tex(i,j-1,1))*anhf(j-1)
     $      +(tex(i-1,j,1)+tex(i,j,1))*ashf(j)) 
     $      + hclu(i,j)*areat(i,j,1)
      enddo
      enddo
c
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
        eni(k)=eni(k)+.5d0*ex(i,j,k)*((u(i,j,k)-sfun(i,j))**2
     $        +(v(i,j,k)-sfvn(i,j))**2)*(areauu(j)*dzu(i,j,k))
        eng =(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k))*ex(i,j,k)
        spm(k)=dmax1(eng,spm(k))
        tmn(k)=tmn(k)+t(i,j,k)*volt(i,j,k)
        smn(k)=smn(k)+s(i,j,k)*volt(i,j,k)
        dmn(k)=dmn(k)+rho(i,j,k)*volt(i,j,k)
        d2mn(k)=d2mn(k)+rho(i,j,k)*rho(i,j,k)*volt(i,j,k)
        eng2=((v(i,j,k)*dzu(i,j,k)+v(i,j-1,k)*dzu(i,j-1,k)
     &      -v(i-1,j,k)*dzu(i-1,j,k)-v(i-1,j-1,k)*dzu(i-1,j-1,k))*dy
     &      -((u(i,j,k)*dzu(i,j,k)+u(i-1,j,k)*dzu(i-1,j,k))*cs(j)
     &      -(u(i,j-1,k)*dzu(i,j-1,k)+u(i-1,j-1,k)*dzu(i-1,j-1,k))
     &      *cs(j-1))*dx)*tex(i,j,k)
        est(k)=est(k)+eng2*eng2*voltr(i,j,k)
!        test1(k) = test1(k) + volt(i,j,k)
!        test2(k) = test2(k) + t(i,j,k)
      enddo
      enddo
      enddo

#ifdef DEBUG1
      do k = 1,km
        write(*,*) nkai,ip,k,tmn(k),smn(k),dmn(k),eni(k),
     &  est(k),spm(k)
!      if(spm(k).eq.0.)then
!         write(*,*)'there is no value -> stop',ip
!         stop
!      endif	  
      enddo
	  write(*,*)
      
!      return
#endif

      call sum_mst(eni,km)
      call sum_mst1(eng,1)
      call sum_mst(tmn,km)
      call sum_mst(smn,km)
      call sum_mst(dmn,km+2)
      call sum_mst(d2mn,km)
      call max_mst(spm,km)
      call max_mst(est,km)

      if(ip.eq.imaster) then
      do k=1,km
        engina=engina+eni(k)/ttvol
        enstro=enstro+.125d0*est(k)/ttvol
        ttmna=ttmna+tmn(k)/ttvol
        ssmna=ssmna+smn(k)/ttvol
        ddmna=ddmna+dmn(k)/ttvol
        eni(k)=eni(k)/tvol(k)
        est(k)=.125d0*est(k)/tvol(k)
        tmn(k)=tmn(k)/tvol(k)
        smn(k)=smn(k)/tvol(k)
        dmn(k)=dmn(k)/tvol(k)
        d2mn(k)=d2mn(k)/tvol(k)
        spm(k)=dsqrt(spm(k))
!       spmx(k)=dsqrt(spmx(k))
      enddo

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

      write(6,*)nkai,'step : mean values'
      write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &     ' hour=',ihour,' min=',imin

      write(*,*)'eng_barotro=',engina,' eng_barocli=',engexa
      write(*,*)'enstrophy=',enstro,' htmean=',htmn,
     &     ' htmean+ice=',htmnpice
      write(*,*) 'level, temp, sal, dens, kin_eng, enstrophy, max_velo'
      do k = 1,km
        write(*,*) k,tmn(k),smn(k),dmn(k),eni(k),est(k),spm(k)
     & ,tvol(k)
      enddo
      write(*,*)'mean t/s/rho = ',ttmna,ssmna,ddmna
      write(*,*)

      endif

      do k = 1,km*7+5
       amax(k) = -100.
       amin(k) = 100.
      enddo
      
      do k = 1,km

      umax = 0.
      vmax = 0.
      tmax = 0.
      smax = 0.
      emax = 0.
      wmax = 0.
      vdtsmax = 0.
      umin = 100.
      vmin = 100.
      tmin = 100.
      smin = 100.
      emin = 100.
      wmin = 100.
      vdtsmin = 100.

      do j = 3,jml-2
      do i = 3,iml-2
      if(umax.lt.u(i,j,k)) then
        umax = u(i,j,k)
        iumax = i
        jumax = j
      endif
      if(umin.gt.u(i,j,k)) then
        umin = u(i,j,k)
        iumin = i
        jumin = j
      endif
      if(vmax.lt.v(i,j,k)) then
        vmax = v(i,j,k)
        ivmax = i
        jvmax = j
      endif
      if(vmin.gt.v(i,j,k)) then
        vmin = v(i,j,k)
        ivmin = i
        jvmin = j
      endif
      if(tex(i,j,k).ne.0.) then
      if(tmax.lt.t(i,j,k)) then
        tmax = t(i,j,k)
        itmax = i
        jtmax = j
      endif
      if(tmin.gt.t(i,j,k)) then
        tmin = t(i,j,k)
        itmin = i
        jtmin = j
      endif
      if(smax.lt.s(i,j,k)) then
        smax = s(i,j,k)
        ismax = i
        jsmax = j
      endif
      if(smin.gt.s(i,j,k)) then
        smin = s(i,j,k)
        ismin = i
        jsmin = j
      endif
      w = wl(i,j,k)/(areat(i,j,k)+1.-tex(i,j,k))
      if(wmax.lt.w) then
        wmax = w
        iwmax = i
        jwmax = j
      endif
      if(wmin.gt.w) then
        wmin = w
        iwmin = i
        jwmin = j
      endif

      if(emax.lt.tke(i,j,k)) then
        emax = tke(i,j,k)
        iemax = i
        jemax = j
      endif
      if(emin.gt.tke(i,j,k)) then
        emin = tke(i,j,k)
        iemin = i
        jemin = j
      endif
      if(vdtsmax.lt.vdts(i,j,k)) then
        vdtsmax = vdts(i,j,k)
        idtmax = i
        jdtmax = j
      endif
      if(vdtsmin.gt.vdts(i,j,k)) then
        vdtsmin = vdts(i,j,k)
        idtmin = i
        jdtmin = j
      endif
      endif

      enddo
      enddo

        amax((k-1)*7+1) = umax
        iamax((k-1)*7+1) = iumax
        jamax((k-1)*7+1) = jumax
        amax((k-1)*7+2) = vmax
        iamax((k-1)*7+2) = ivmax
        jamax((k-1)*7+2) = jvmax
        amax((k-1)*7+3) = tmax
        iamax((k-1)*7+3) = itmax
        jamax((k-1)*7+3) = jtmax
        amax((k-1)*7+4) = smax
        iamax((k-1)*7+4) = ismax
        jamax((k-1)*7+4) = jsmax
        amax((k-1)*7+5) = wmax
        iamax((k-1)*7+5) = iwmax
        jamax((k-1)*7+5) = jwmax
        amax((k-1)*7+6) = emax
        iamax((k-1)*7+6) = iemax
        jamax((k-1)*7+6) = jemax
        amax((k-1)*7+7) = vdtsmax
        iamax((k-1)*7+7) = idtmax
        jamax((k-1)*7+7) = jdtmax
c
        amin((k-1)*7+1) = umin
        iamin((k-1)*7+1) = iumin
        jamin((k-1)*7+1) = jumin
        amin((k-1)*7+2) = vmin
        iamin((k-1)*7+2) = ivmin
        jamin((k-1)*7+2) = jvmin
        amin((k-1)*7+3) = tmin
        iamin((k-1)*7+3) = itmin
        jamin((k-1)*7+3) = jtmin
        amin((k-1)*7+4) = smin
        iamin((k-1)*7+4) = ismin
        jamin((k-1)*7+4) = jsmin
        amin((k-1)*7+5) = wmin
        iamin((k-1)*7+5) = iwmin
        jamin((k-1)*7+5) = jwmin
        amin((k-1)*7+6) = emin
        iamin((k-1)*7+6) = iemin
        jamin((k-1)*7+6) = jemin
        amin((k-1)*7+7) = vdtsmin
        iamin((k-1)*7+7) = idtmin
        jamin((k-1)*7+7) = jdtmin
c
      enddo
c
      aimax = 0.
      himax = 0.
      volimax=0.
      uimax = 0.
      vimax = 0.
      uimin = 0.
      vimin = 0.
      hclmax=-1.d3
      hclmin= 1.d3
      ummax=0.d0
      vmmax=0.d0
      ummin=0.d0
      vmmin=0.d0
c
      do j=3,jml-2
      do i=3,iml-2
      if(aice(i,j).gt.0.d0)then
         if(aimax.lt.aice(i,j)) then
            aimax = aice(i,j)
            iaimax = i
            jaimax = j
         endif
         if(volimax.lt.volice(i,j)) then
            volimax = volice(i,j)
            ivolimax = i
            jvolimax = j
         endif
         if(himax.lt.volice(i,j)/areat(i,j,1)/aice(i,j)) then
            himax = volice(i,j)/areat(i,j,1)/aice(i,j)
            ihimax = i
            jhimax = j
         endif
         if(uimax.lt.uice(i,j)) then
            uimax = uice(i,j)
            iuimax = i
            juimax = j
         endif
         if(uimin.gt.uice(i,j)) then
            uimin = uice(i,j)
            iuimin = i
            juimin = j
         endif
          if(vimax.lt.vice(i,j)) then
            vimax = vice(i,j)
            ivimax = i
            jvimax = j
         endif
         if(vimin.gt.vice(i,j)) then
            vimin = vice(i,j)
            ivimin = i
            jvimin = j
         endif
      endif
c
       if(tex(i,j,1).eq.1.) then
         if(hclmax.lt.hcl(i,j)) then
            hclmax = hcl(i,j)
            ihclmax = i
            jhclmax = j
         endif
         if(hclmin.gt.hcl(i,j)) then
            hclmin = hcl(i,j)
            ihclmin = i
            jhclmin = j
         endif
       endif
         if(ummax.lt.sfun(i,j)) then
            ummax = sfun(i,j)
            iummax = i
            jummax = j
         endif
         if(ummin.gt.sfun(i,j)) then
            ummin = sfun(i,j)
            iummin = i
            jummin = j
         endif
          if(vmmax.lt.sfvn(i,j)) then
            vmmax = sfvn(i,j)
            ivmmax = i
            jvmmax = j
         endif
         if(vmmin.gt.sfvn(i,j)) then
            vmmin = sfvn(i,j)
            ivmmin = i
            jvmmin = j
         endif
      enddo
      enddo
c
        amax(km*7+1) = aimax
        iamax(km*7+1) = iaimax
        jamax(km*7+1) = jaimax
        amax(km*7+2) = himax
        iamax(km*7+2) = ihimax
        jamax(km*7+2) = jhimax
        amax(km*7+3) = volimax
        iamax(km*7+3) = ivolimax
        jamax(km*7+3) = jvolimax
        amax(km*7+4) = uimax
        iamax(km*7+4) = iuimax
        jamax(km*7+4) = juimax
        amax(km*7+5) = vimax
        iamax(km*7+5) = ivimax
        jamax(km*7+5) = jvimax
        amin(km*7+4) = uimin
        iamin(km*7+4) = iuimin
        jamin(km*7+4) = juimin
        amin(km*7+5) = vimin
        iamin(km*7+5) = ivimin
        jamin(km*7+5) = jvimin
        amax(km*7+6) = hclmax
        iamax(km*7+6) = ihclmax
        jamax(km*7+6) = jhclmax
        amax(km*7+7) = ummax
        iamax(km*7+7) = iummax
        jamax(km*7+7) = jummax
        amax(km*7+8) = vmmax
        iamax(km*7+8) = ivmmax
        jamax(km*7+8) = jvmmax
        amin(km*7+6) = hclmin
        iamin(km*7+6) = ihclmin
        jamin(km*7+6) = jhclmin
        amin(km*7+7) = ummin
        iamin(km*7+7) = iummin
        jamin(km*7+7) = jummin
        amin(km*7+8) = vmmin
        iamin(km*7+8) = ivmmin
        jamin(km*7+8) = jvmmin
c
      call max_at_2d(amax,iamax,jamax,km*7+8)
      call min_at_2d(amin,iamin,jamin,km*7+8)
c
      if(ip.eq.imaster) then
      write(*,*) nkai,'step : max and min values'
c
#ifdef EQ100M
      do k = 1,5
#else
      do k = 1,km
#endif
      umax  = amax((k-1)*7+1)
      iumax = iamax((k-1)*7+1)
      jumax = jamax((k-1)*7+1)
      vmax  = amax((k-1)*7+2)
      ivmax = iamax((k-1)*7+2)
      jvmax = jamax((k-1)*7+2)
      tmax  = amax((k-1)*7+3)
      itmax = iamax((k-1)*7+3)
      jtmax = jamax((k-1)*7+3)
      smax  = amax((k-1)*7+4)
      ismax = iamax((k-1)*7+4)
      jsmax = jamax((k-1)*7+4)
      wmax  = amax((k-1)*7+5)
      iwmax = iamax((k-1)*7+5)
      jwmax = jamax((k-1)*7+5)
      emax  = amax((k-1)*7+6)
      iemax = iamax((k-1)*7+6)
      jemax = jamax((k-1)*7+6)
      vdtsmax  = amax((k-1)*7+7)
      idtmax = iamax((k-1)*7+7)
      jdtmax = jamax((k-1)*7+7)
c
      umin  = amin((k-1)*7+1)
      iumin = iamin((k-1)*7+1)
      jumin = jamin((k-1)*7+1)
      vmin  = amin((k-1)*7+2)
      ivmin = iamin((k-1)*7+2)
      jvmin = jamin((k-1)*7+2)
      tmin  = amin((k-1)*7+3)
      itmin = iamin((k-1)*7+3)
      jtmin = jamin((k-1)*7+3)
      smin  = amin((k-1)*7+4)
      ismin = iamin((k-1)*7+4)
      jsmin = jamin((k-1)*7+4)
      wmin  = amin((k-1)*7+5)
      iwmin = iamin((k-1)*7+5)
      jwmin = jamin((k-1)*7+5)
      emin  = amin((k-1)*7+6)
      iemin = iamin((k-1)*7+6)
      jemin = jamin((k-1)*7+6)
      vdtsmin  = amin((k-1)*7+7)
      idtmin = iamin((k-1)*7+7)
      jdtmin = jamin((k-1)*7+7)
c
      write(*,*) 'k = ',k,
     &     ' (dp= ',nint(dp(k)),' dz= ',nint(dz(k)),'cm)'
      write(*,701)iumax,jumax,umax,iumin,jumin,umin,
     &     ivmax,jvmax,vmax,ivmin,jvmin,vmin
      write(*,702)itmax,jtmax,tmax,itmin,jtmin,tmin,
     &     ismax,jsmax,smax,ismin,jsmin,smin
      write(*,703)iwmax,jwmax,wmax,iwmin,jwmin,wmin
      write(*,704)iemax,jemax,emax,iemin,jemin,emin,
     &     idtmax,jdtmax,vdtsmax,idtmin,jdtmin,vdtsmin
c
      enddo
c
      write(*,*)
c
        aimax = amax(km*7+1)
        iaimax  = iamax(km*7+1)
        jaimax = jamax(km*7+1)
        himax = amax(km*7+2)
        ihimax = iamax(km*7+2)
        jhimax = jamax(km*7+2)
        volimax = amax(km*7+3)
        ivolimax = iamax(km*7+3)
        jvolimax = jamax(km*7+3)
        uimax = amax(km*7+4)
        iuimax = iamax(km*7+4)
        juimax = jamax(km*7+4)
        vimax = amax(km*7+5)
        ivimax = iamax(km*7+5)
        jvimax = jamax(km*7+5)
        uimin = amin(km*7+4)
        iuimin = iamin(km*7+4)
        juimin = jamin(km*7+4)
        vimin = amin(km*7+5)
        ivimin = iamin(km*7+5)
        jvimin = jamin(km*7+5)
c
      write(*,*)nkai,'step : ice'
      write(*,705)iaimax,jaimax,aimax,
     &     ihimax,jhimax,himax,ivolimax,jvolimax,volimax
      write(*,706)iuimax,juimax,uimax,iuimin,juimin,uimin,
     &     ivimax,jvimax,vimax,ivimin,jvimin,vimin
      write(*,*)
      endif
c
#ifdef PC68M
 701  format(' u: max(',i4,' ',i4,')',f7.2,', min(',i4,' ',i4,')',f7.2,
     &     '   v: max(',i4,' ',i4,')',f7.2,', min(',i4,' ',i4,')',f7.2)
 702  format(' t: max(',i4,' ',i4,')',f7.3,', min(',i4,' ',i4,')',f7.3,
     &     '   s: max(',i4,' ',i4,')',f7.3,', min(',i4,' ',i4,')',f7.3)
 703  format(' w: max(',i4,' ',i4,')',f10.6,
     &     ', min(',i4,' ',i4,')',f10.6)
 704  format(' tke: max(',i4,' ',i4,')',g10.4,
     &     ', min(',i4,' ',i4,')',g10.4,
     &     '   vdts: max(',i4,' ',i4,')',f8.3,
     &     ', min(',i4,' ',i4,')',f8.3)
 705  format(' aice: max(',i4,' ',i4,')',f7.4,
     &     '   hice: max(',i4,' ',i4,')',f7.2,
     &     '   volice: max(',i4,' ',i4,')',g15.4)
 706  format(' uice: max(',i4,' ',i4,')',f7.2,
     &     ', min(',i4,' ',i4,')',f7.2,
     &     '   vice: max(',i4,' ',i4,')',f7.2,
     &     ', min(',i4,' ',i4,')',f7.2)
#else
 701  format(' u: max(',i3,' ',i3,')',f7.2,', min(',i3,' ',i3,')',f7.2,
     &     '   v: max(',i3,' ',i3,')',f7.2,', min(',i3,' ',i3,')',f7.2)
 702  format(' t: max(',i3,' ',i3,')',f7.3,', min(',i3,' ',i3,')',f7.3,
     &     '   s: max(',i3,' ',i3,')',f7.3,', min(',i3,' ',i3,')',f7.3)
 703  format(' w: max(',i3,' ',i3,')',f10.6,
     &     ', min(',i3,' ',i3,')',f10.6)
 704  format(' tke: max(',i3,' ',i3,')',g10.4,
     &     ', min(',i3,' ',i3,')',g10.4,
     &     '   vdts: max(',i3,' ',i3,')',f8.3,
     &     ', min(',i3,' ',i3,')',f8.3)
 705  format(' aice: max(',i3,' ',i3,')',f7.4,
     &     '   hice: max(',i3,' ',i3,')',f7.2,
     &     '   volice: max(',i3,' ',i3,')',g15.4)
 706  format(' uice: max(',i3,' ',i3,')',f7.2,
     &     ', min(',i3,' ',i3,')',f7.2,
     &     '   vice: max(',i3,' ',i3,')',f7.2,
     &     ', min(',i3,' ',i3,')',f7.2)
#endif
c
      hclmax  =  amax(km*7+6)
      ihclmax = iamax(km*7+6)
      jhclmax = jamax(km*7+6)
      ummax   =  amax(km*7+7)
      iummax  = iamax(km*7+7)
      jummax  = jamax(km*7+7)
      vmmax   =  amax(km*7+8)
      ivmmax  = iamax(km*7+8)
      jvmmax  = jamax(km*7+8)
      hclmin  =  amin(km*7+6)
      ihclmin = iamin(km*7+6)
      jhclmin = jamin(km*7+6)
      ummin   =  amin(km*7+7)
      iummin  = iamin(km*7+7)
      jummin  = jamin(km*7+7)
      vmmin   =  amin(km*7+8)
      ivmmin  = iamin(km*7+8)
      jvmmin  = jamin(km*7+8)
c
      if(ip.eq.imaster) then
      write(*,*)nkai,'step : brtrop'
      write(*,707) ihclmax,jhclmax,hclmax,ihclmin,jhclmin,hclmin
      write(*,708) iummax,jummax,ummax,iummin,jummin,ummin,
     &  ivmmax,jvmmax,vmmax,ivmmin,jvmmin,vmmin
      write(*,*)
c
      endif
c
#ifdef PC68M
 707  format(' hcl: max(',i4,' ',i4,')',f8.2,
     &     ', min(',i4,' ',i4,')',f8.2)
 708  format(' sfu: max(',i4,' ',i4,')',f7.2,
     &     ', min(',i4,' ',i4,')',f7.2),
     &     ' sfv: max(',i4,' ',i4,')',f7.2,
     &     ', min(',i4,' ',i4,')',f7.2)
#else
 707  format(' hcl: max(',i3,' ',i3,')',f7.2,
     &     ', min(',i3,' ',i3,')',f7.2)
 708  format(' sfu: max(',i3,' ',i3,')',f7.2,
     &     ', min(',i3,' ',i3,')',f7.2),
     &     ' sfv: max(',i3,' ',i3,')',f7.2,
     &     ', min(',i3,' ',i3,')',f7.2)
#endif
c
c
      if(ip.eq.imaster.and.spm(1).eq.0.)then
         write(*,*)'there is no value -> stop'
         stop
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine spinup                                             c
c
c     toyoda
c                                                                   c
c     Spin-up diagnosation.                                         c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spinup
      use param
      use mod_mpi
      implicit none
#include "common.h"
	  
      real*8 engall,tr_kuro,tr_oya,tr_ekm,tr_tsugaru,tr_wtsm,tr_etsm
      real*8 pht24n100,pht24n1000,pht24nall
      real*8 fcx,fcy,s1,s2
      real*8 totalivol,totaliarea
      integer i,j,k,l,iia,jja1,jja2,jja,iia1,iia2

      engexa=0.d0
      do j=3,jml-2
      do i=3,iml-2
        engexa=engexa+.5d0*(sfun(i,j)*sfun(i,j)+sfvn(i,j)*sfvn(i,j))
     1    *areauu(j)*(hrr(i,j)+hclu(i,j))
      enddo
      enddo
      call sum_mst1(engexa,1)
      if(ip.eq.imaster) then
         engexa=engexa/ttvol
      endif
c
      do k=1,km
         eni(k)=0.
      enddo

      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
        eni(k)=eni(k)+.5d0*ex(i,j,k)*((u(i,j,k)-sfun(i,j))**2
     $        +(v(i,j,k)-sfvn(i,j))**2)*(areauu(j)*dzu(i,j,k))
      enddo
      enddo
      enddo
      call sum_mst(eni,km)
c
      engina=0.
      if(ip.eq.imaster) then
      do k=1,km
        engina=engina+eni(k)/ttvol
      enddo
      endif
c
      engall=0.
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
        engall=engall+.5d0*ex(i,j,k)
     &        *( u(i,j,k)**2+v(i,j,k)**2 )*areauu(j)*dzu(i,j,k)
      enddo
      enddo
      enddo
      call sum_mst1(engall,1)
      engall=engall/ttvol

      tr_kuro=0.
      fcx=2./sqrt(5.)
      fcy=1./sqrt(5.)
      do l=1,30
         iia=161+(l-1)
         jja1=353-2*(l-1)
         jja2=353-2*(l-1)-1

         do j=3,jml-2
         do i=3,iml-2
            if(lw(ip_x)+i-1.eq.iia .and. ls(ip_y)+j-1.eq.jja1)then
              s1=sqrt((dx*cs(j))**2+dy**2)
            do k=1,67 !1000m
               tr_kuro=tr_kuro
     &           +dmax1(fcx*u(iia,jja1,k)+fcy*v(iia,jja1,k),0.d0)
     &           *ex(iia,jja1,k)*s1*dzu(iia,jja1,k)
     &           +dmax1(fcx*u(iia,jja2,k)+fcy*v(iia,jja2,k),0.d0)
     &           *ex(iia,jja2,k)*s1*dzu(iia,jja2,k)
c            write(*,*)l,iia,jja1,jja2,k,s1,s2,tr_kuro
            enddo
            endif
         enddo
         enddo
      enddo
      call sum_mst1(tr_kuro,1)

      tr_oya=0.
      jja=426
      iia1=234
      iia2=264
      do j=3,jml-2
      do i=3,iml-2
         if(lw(ip_x)+i-1.ge.iia1 .and. lw(ip_x)+i-1.le.iia2) then
         if(ls(ip_y)+j-1.eq.jja)then
         do k=1,67
            tr_oya=tr_oya
     &           +dmin1(v(i,j,k),0.d0)*ex(i,j,k)*dx*cs(j)*dzu(i,j,k)
         enddo
         endif
         endif
      enddo
      enddo
      call sum_mst1(tr_oya,1)

      tr_ekm=0.
      jja=522
      iia1=340
      iia2=370
      do j=3,jml-2
      do i=3,iml-2
         if(lw(ip_x)+i-1.ge.iia1 .and. lw(ip_x)+i-1.le.iia2) then
         if(ls(ip_y)+j-1.eq.jja)then
         do k=1,67
            tr_ekm=tr_ekm
     &           +dmin1(v(i,j,k),0.d0)*ex(i,j,k)*dx*cs(j)*dzu(i,j,k)
         enddo
         endif
         endif
      enddo
      enddo
      call sum_mst1(tr_ekm,1)

      pht24n100=0.
      jja=274
      do j=3,jml-2
         if(ls(ip_y)+j-1.eq.jja)then
         do k=1,21
         do i=3,iml-2
            pht24n100=pht24n100
     &           +v(i,j,k)*ex(i,j,k)*dx*cs(j)*dzu(i,j,k)
     &           *(t(i,j,k)+t(i+1,j,k)+t(i,j+1,k)+t(i+1,j+1,k))/4.
         enddo
         enddo
         endif
      enddo
      call sum_mst1(pht24n100,1)

      pht24n1000=0.
      jja=274
      do j=3,jml-2
         if(ls(ip)+j-1.eq.jja)then
         do k=1,67
         do i=3,im-2
            pht24n1000=pht24n1000
     &           +v(i,j,k)*ex(i,j,k)*dx*cs(j)*dzu(i,j,k)
     &           *(t(i,j,k)+t(i+1,j,k)+t(i,j+1,k)+t(i+1,j+1,k))/4.
         enddo
         enddo
         endif
      enddo
      call sum_mst1(pht24n1000,1)

      pht24nall=0.
      jja=274
      do j=3,jml-2
         if(ls(ip)+j-1.eq.jja)then
         do k=1,km
         do i=3,iml-2
            pht24nall=pht24nall
     &           +v(i,j,k)*ex(i,j,k)*dx*cs(j)*dzu(i,j,k)
     &           *(t(i,j,k)+t(i+1,j,k)+t(i,j+1,k)+t(i+1,j+1,k))/4.
         enddo
         enddo
         endif
      enddo
      call sum_mst1(pht24nall,1)

      totalivol=0.d0
      do i=3,iml-2
      do j=3,jml-2
         totalivol=totalivol+areat(i,j,1)*volice(i,j)
      enddo
      enddo
      call sum_mst1(totalivol,1)
c
      totaliarea=0.d0
      do i=3,iml-2
      do j=3,jml-2
         totaliarea=totaliarea+areat(i,j,1)*aice(i,j)
      enddo
      enddo
      call sum_mst1(totaliarea,1)

      tr_tsugaru=0.
      jja1=413
      jja2=417
      iia=216

      do j=3,jml-2
      do i=3,iml-2
         if(ls(ip_y)+j-1.ge.jja1 .and. ls(ip_y)+j-1.le.jja2) then
         if(lw(ip_x)+i-1.eq.iia)then
         do k=1,km
            tr_tsugaru=tr_tsugaru
     &           +u(i,j,k)*ex(i,j,k)*dx*cs(j)*dzu(i,j,k)
         enddo
         endif
         endif
      enddo
      enddo
      call sum_mst1(tr_tsugaru,1)

      fcx=1./sqrt(5.)
      fcy=2./sqrt(5.)
      tr_wtsm=0.
      do l=1,3
         iia1=146-2*(l-1)
         iia2=146-2*(l-1)-1
         jja=357+(l-1)

         do j=3,jml-2
         do i=3,iml-2
            if(lw(ip_x)+i-1.eq.iia1 .and. ls(ip_y)+j-1.eq.jja)then
              s1=sqrt((dx*cs(j))**2+dy**2)
            do k=1,km
               tr_wtsm=tr_wtsm
     &           +(fcx*u(iia1,jja,k)+fcy*v(iia1,jja,k))
     &           *ex(iia1,jja,k)*s1*dzu(iia1,jja,k)
     &           +(fcx*u(iia2,jja,k)+fcy*v(iia2,jja,k))
     &           *ex(iia2,jja,k)*s1*dzu(iia2,jja,k)
            enddo
            endif
         enddo
         enddo
      enddo
      call sum_mst1(tr_wtsm,1)
      
      tr_etsm=0.
      do l=1,5
         iia1=155-2*(l-1)
         iia2=155-2*(l-1)-1
         jja=352+(l-1)

         do j=3,jml-2
         do i=3,iml-2
            if(lw(ip_x)+i-1.eq.iia1 .and. ls(ip_y)+j-1.eq.jja)then
              s1=sqrt((dx*cs(j))**2+dy**2)
            do k=1,km
               tr_etsm=tr_etsm
     &           +(fcx*u(iia1,jja,k)+fcy*v(iia1,jja,k))
     &           *ex(iia1,jja,k)*s1*dzu(iia1,jja,k)
     &           +(fcx*u(iia2,jja,k)+fcy*v(iia2,jja,k))
     &           *ex(iia2,jja,k)*s1*dzu(iia2,jja,k)
            enddo
            endif
         enddo
         enddo
      enddo
      call sum_mst1(tr_etsm,1)

      if(ip.eq.imaster)then
         write(spinupf)ahour,engexa,engina,engall,
     &        tr_kuro,tr_oya,tr_ekm,pht24n100,pht24n1000,pht24nall,
     &        totalivol,totaliarea,tr_tsugaru,tr_wtsm,tr_etsm
      endif
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine vort_eg                                            c
c
c     toyoda
c                                                                   c
c     Term ballance of vorticity equation in the EG region          c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef CHECKVORT
      subroutine vort_eg

      use param
      use mod_mpi
      implicit none
#include "common.h"
c
      real*8 zeta(im,jm,0:km+1)
      real*8 circ,adv_w,adv_e,adv_s,adv_n,curltau,betav,btmstr
      real*8 zdwdz,fdwdz,dwdxdvdz,dwdydudz,allterms,adv
      real*8 visc_lp,visc_bh
      real*8 faca,rho0,beta(jm),vel00,vel10,vel01,vel11,hduv_bh_x
      real*8 zeta2x(im,jm,0:km+1),zeta2y(im,jm,0:km+1)
c      real*8 visc_lp2,uafr00,uafr10,uafr01,uafr11,
c     &     vafr00,vafr10,vafr01,vafr11
      real*8 betav2w,betav2e,betav2s,betav2n,betav2,facv
c
      if(mod(nkai,nxday).eq.0.and.matsno.eq.0)then
c
c
      ii1=613
      ii2=691
      jj1=219
      jj2=283
c      ii1=613
c      ii2=614
c      jj1=219
c      jj2=220
c
      rho0=1.03
      do j=2,jml-1
         beta(j)=1./radius*2.*omega*cs(j)
      enddo
      hduv_bh_x=1.d17
c
c
      call exch_t3d_w1p(u,v,1,2)
      call exch_t3d_s1p(u,v,3,4)
      call exch_t2d_w1p(wsx,wsy,9)
      call exch_t2d_s1p(wsx,wsy,10)
c
      do i=1,im
      do j=1,jm
      do k=0,km+1
         zeta(i,j,k)=0.
         zeta2x(i,j,k)=0.
         zeta2y(i,j,k)=0.
      enddo
      enddo
      enddo
c
      call wait_t3d_w1p(u,v,1,2)
      call wait_t3d_s1p(u,v,3,4)
c
      do i=3,iml-2
      do j=3,jml-2
      do k=1,km
c         zeta(i,j,k)=tex(i,j,k)*(
c     &        (v(i,j-1,k)-v(i-1,j-1,k))/dx/cs(j-1)/2.
c     &        +(v(i,j,k)-v(i-1,j,k))/dx/cs(j)/2.
c     &        -(u(i-1,j,k)-u(i-1,j-1,k))/dy/2.
c     &        -(u(i,j,k)-u(i,j-1,k))/dy/2.
c     &        )
         zeta(i,j,k)=tex(i,j,k)*(
     &        ((v(i,j,k)+v(i,j-1,k))/(1.+ex(i,j,k)*ex(i,j-1,k))
     &        -(v(i-1,j,k)+v(i-1,j-1,k))/(1.+ex(i-1,j,k)*ex(i-1,j-1,k)))
     &        /dx/cst(j)
     &        -((u(i,j,k)+u(i-1,j,k))/(1.+ex(i,j,k)*ex(i-1,j,k))
     &        -(u(i,j-1,k)+u(i-1,j-1,k))/(1.+ex(i,j-1,k)*ex(i-1,j-1,k)))
     &        /dy )
c         if(lw(ip_x)+i-1.eq.(ii1+112)/2.and.
c     &        ls(ip_y)+j-1.eq.(jj1+jj2)/2.)
c     &        write(*,*)lw(ip_x)+i-1,ls(ip_y)+j-1,k,zeta(i,j,k)
      enddo
      enddo
      enddo
c      stop
c
      call exch_t3d_w1(zeta,5)
      call exch_t3d_e1(zeta,6)
      call exch_t3d_s1(zeta,7)
      call exch_t3d_n1(zeta,8)
c
      circ=0.
      do i=3,iml-2
      do j=3,jml-2
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2.and.
     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2)then
c         faca=(anhf(j-1)+ashf(j))*2.
c         if(ls(ip_y)+j-1.eq.jj1)faca=ashf(j)*2.
c         if(ls(ip_y)+j-1.eq.jj2)faca=anhf(j-1)*2.
c         if(lw(ip_x)+i-1.eq.ii1.or.lw(ip_x)+i-1.eq.ii2)faca=faca/2.
         do k=1,km
c            circ=circ+zeta(i,j,k)*dzt(i,j,k)*faca
c
            facv=0.
            if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
     &           facv=facv+ex(i-1,j-1,k)*anhf(j-1)*dzu(i-1,j-1,k)
            if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
     &           facv=facv+ex(i-1,j,k)*ashf(j)*dzu(i-1,j,k)
            if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
     &           facv=facv+ex(i,j-1,k)*anhf(j-1)*dzu(i,j-1,k)
            if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
     &           facv=facv+ex(i,j,k)*ashf(j)*dzu(i,j,k)
            circ=circ+zeta(i,j,k)*facv
c
c         if(lw(ip_x)+i-1.eq.(ii1+112)/2.and.
c     &        ls(ip_y)+j-1.eq.(jj1+jj2)/2.)
c     &        write(*,*)lw(ip_x)+i-1,ls(ip_y)+j-1,k,zeta(i,j,k)*facv
c
         enddo
      endif
      enddo
      enddo
      call sum_mst1(circ,1)
c
      adv_w=0.
      do i=3,iml-2
      do j=3,jml-2
      do k=1,km
      if(lw(ip_x)+i-1.eq.ii1.and.
     &     ls(ip_y)+j-1.ge.jj1+1.and.ls(ip_y)+j-1.le.jj2)
     &     adv_w=adv_w+zeta(i,j,k)*(u(i-1,j-1,k)+u(i,j-1,k))/2.
     &     *dzumin(i,j-1,k)*dy/2.
      if(lw(ip_x)+i-1.eq.ii1.and.
     &     ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2-1)
     &     adv_w=adv_w+zeta(i,j,k)*(u(i-1,j,k)+u(i,j,k))/2.
     &     *dzumin(i,j,k)*dy/2.
      enddo
      enddo
      enddo
      call sum_mst1(adv_w,1)
c
      adv_e=0.
      do i=3,iml-2
      do j=3,jml-2
      do k=1,km
      if(lw(ip_x)+i-1.eq.ii2.and.
     &     ls(ip_y)+j-1.ge.jj1+1.and.ls(ip_y)+j-1.le.jj2)
     &     adv_e=adv_e-zeta(i,j,k)*(u(i-1,j-1,k)+u(i,j-1,k))/2.
     &     *dzumin(i,j-1,k)*dy/2.
      if(lw(ip_x)+i-1.eq.ii2.and.
     &     ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2-1)
     &     adv_e=adv_e-zeta(i,j,k)*(u(i-1,j,k)+u(i,j,k))/2.
     &     *dzumin(i,j,k)*dy/2.
      enddo
      enddo
      enddo
      call sum_mst1(adv_e,1)
c
      adv_s=0.
      do i=3,iml-2
      do j=3,jml-2
      do k=1,km
      if(lw(ip_x)+i-1.ge.ii1+1.and.lw(ip_x)+i-1.le.ii2
     &        .and.ls(ip_y)+j-1.eq.jj1)
     &        adv_s=adv_s+zeta(i,j,k)*(v(i-1,j-1,k)+v(i-1,j,k))/2.
     &        *dzvmin(i-1,j,k)*dx*cst(j)/2.
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2-1
     &        .and.ls(ip_y)+j-1.eq.jj1)
     &        adv_s=adv_s+zeta(i,j,k)*(v(i,j-1,k)+v(i,j,k))/2.
     &        *dzvmin(i,j,k)*dx*cst(j)/2.
      enddo
      enddo
      enddo
      call sum_mst1(adv_s,1)
c
      adv_n=0.
      do i=3,iml-2
      do j=3,jml-2
      do k=1,km
      if(lw(ip_x)+i-1.ge.ii1+1.and.lw(ip_x)+i-1.le.ii2
     &        .and.ls(ip_y)+j-1.eq.jj2)
     &        adv_n=adv_n-zeta(i,j,k)*(v(i-1,j-1,k)+v(i-1,j,k))/2.
     &        *dzvmin(i-1,j,k)*dx*cst(j)/2.
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2-1
     &        .and.ls(ip_y)+j-1.eq.jj2)
     &        adv_n=adv_n-zeta(i,j,k)*(v(i,j-1,k)+v(i,j,k))/2.
     &        *dzvmin(i,j,k)*dx*cst(j)/2.
      enddo
      enddo
      enddo
      call sum_mst1(adv_n,1)
c
      call wait_t2d_w1p(wsx,wsy,9)
      call wait_t2d_s1p(wsx,wsy,10)
c
      curltau=0.
      do i=3,iml-2
      do j=3,jml-2
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2.and.
     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2)then
         faca=0.
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        faca=faca+ex(i-1,j-1,1)*anhf(j-1)
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
     &        faca=faca+ex(i-1,j,1)*ashf(j)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        faca=faca+ex(i,j-1,1)*anhf(j-1)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
     &        faca=faca+ex(i,j,1)*ashf(j)
c         curltau=curltau
c     &        +( (wsy(i,j)-wsy(i-1,j))/dx/cs(j)/2.
c     &        +(wsy(i,j-1)-wsy(i-1,j-1))/dx/cs(j-1)/2.
c     &        -(wsx(i,j)-wsx(i,j-1))/dy/2.
c     &        -(wsx(i-1,j)-wsx(i-1,j-1))/dy/2. )
c     &        /rho0*faca
         curltau=curltau+tex(i,j,1)*(
     &        +((wsy(i,j)+wsy(i,j-1))/(1.+ex(i,j,1)*ex(i,j-1,1))
     &        -(wsy(i-1,j)+wsy(i-1,j-1))/(1.+ex(i-1,j,1)*ex(i-1,j-1,1)))
     &        /dx/cst(j)
     &        -((wsx(i,j)+wsx(i-1,j))/(1.+ex(i,j,1)*ex(i-1,j,1))
     &        -(wsx(i,j-1)+wsx(i-1,j-1))/(1.+ex(i,j-1,1)*ex(i-1,j-1,1)))
     &        /dy
     &        )/rho0*faca
      endif
      enddo
      enddo
      call sum_mst1(curltau,1)
c
c      betav=0.
      betav2w=0.
      betav2e=0.
      betav2s=0.
      betav2n=0.
      do i=3,iml-2
      do j=3,jml-2
c      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2-1.and.
c     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2-1)then
c         do k=1,km
c            betav=betav
c     &           -v(i,j,k)*beta(j)*dzu(i,j,k)*areauu(j)*ex(i,j,k)
c         enddo
c      endif
      do k=1,km
      if(lw(ip_x)+i-1.eq.ii1.and.
     &     ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2-1)
     &     betav2w=betav2w
     &     +(u(i-1,j,k)+u(i,j,k))/2.*cor(j)*dzumin(i,j,k)*dy
      if(lw(ip_x)+i-1.eq.ii2.and.
     &     ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2-1)
     &     betav2e=betav2e
     &     -(u(i-1,j,k)+u(i,j,k))/2.*cor(j)*dzumin(i,j,k)*dy
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2-1.and.
     &     ls(ip_y)+j-1.eq.jj1)
     &     betav2s=betav2s
     &     +(v(i,j-1,k)+v(i,j,k))/2.*(cor(j-1)+cor(j))/2.
     &     *dzvmin(i,j,k)*dx*cst(j)
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2-1.and.
     &     ls(ip_y)+j-1.eq.jj2)
     &     betav2n=betav2n
     &     -(v(i,j-1,k)+v(i,j,k))/2.*(cor(j-1)+cor(j))/2.
     &     *dzvmin(i,j,k)*dx*cst(j)
      enddo
      enddo
      enddo
c      call sum_mst(betav,1)
      call sum_mst1(betav2w,1)
      call sum_mst1(betav2e,1)
      call sum_mst1(betav2s,1)
      call sum_mst1(betav2n,1)
c
      call wait_t3d_w1(zeta,5)
      call wait_t3d_e1(zeta,6)
      call wait_t3d_s1(zeta,7)
      call wait_t3d_n1(zeta,8)
c
      do i=3,iml-2
      do j=3,jml-2
      do k=1,km
c         zeta2x(i,j,k)=(zeta(i+1,j,k)-2.*zeta(i,j,k)+zeta(i-1,j,k))
c     &           /dx/cst(j)/dx/cst(j)
c         zeta2y(i,j,k)=(zeta(i,j+1,k)-2.*zeta(i,j,k)+zeta(i,j-1,k))
c     &           /dy/dy
         zeta2x(i,j,k)=tex(i,j,k)*
     &        ( (zeta(i+1,j,k)-zeta(i,j,k))*tex(i+1,j,k)
     &        -(zeta(i,j,k)-zeta(i-1,j,k))*tex(i-1,j,k) )
     &        /dx/cst(j)/dx/cst(j)
         zeta2y(i,j,k)=tex(i,j,k)*
     &        ( (zeta(i,j+1,k)-zeta(i,j,k))*tex(i,j+1,k)
     &        -(zeta(i,j,k)-zeta(i,j-1,k))*tex(i,j-1,k) )
     &        /dy/dy
      enddo
      enddo
      enddo
      call exch_t3d_w1(zeta2x,11)
      call exch_t3d_e1(zeta2x,12)
      call exch_t3d_s1(zeta2y,13)
      call exch_t3d_n1(zeta2y,14)
c
      visc_lp=0.
c      visc_lp2=0.
      visc_lp3=0.
      do i=3,iml-2
      do j=3,jml-2
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2.and.
     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2)then
c         faca=0.
c         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
c     &        faca=faca+ex(i-1,j-1,1)*anhf(j-1)
c         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
c     &        faca=faca+ex(i-1,j,1)*ashf(j)
c         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
c     &        faca=faca+ex(i,j-1,1)*anhf(j-1)
c         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
c     &        faca=faca+ex(i,j,1)*ashf(j)
         do k=1,km
c            visc_lp=visc_lp+hduv*(
c     &           (zeta(i+1,j,k)-2.*zeta(i,j,k)+zeta(i-1,j,k))
c     &           /dx/cst(j)/dx/cst(j)
c     &           +(zeta(i,j+1,k)-2.*zeta(i,j,k)+zeta(i,j-1,k))
c     &           /dy/dy )*faca*dzt(i,j,k)
c
         facv=0.
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        facv=facv+ex(i-1,j-1,k)*anhf(j-1)*dzu(i-1,j-1,k)
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
     &        facv=facv+ex(i-1,j,k)*ashf(j)*dzu(i-1,j,k)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        facv=facv+ex(i,j-1,k)*anhf(j-1)*dzu(i,j-1,k)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
     &        facv=facv+ex(i,j,k)*ashf(j)*dzu(i,j,k)
            visc_lp=visc_lp+tex(i,j,k)*hduv*(
     &           ( (zeta(i+1,j,k)-zeta(i,j,k))*tex(i+1,j,k)
     &           -(zeta(i,j,k)-zeta(i-1,j,k))*tex(i-1,j,k) )
     &           /dx/cst(j)/dx/cst(j)
     &           +( (zeta(i,j+1,k)-zeta(i,j,k))*tex(i,j+1,k)
     &           -(zeta(i,j,k)-zeta(i,j-1,k))*tex(i,j-1,k) )
     &           /dy/dy )*facv
c
c            uafr00=(u(i,j-1,k)-2.*u(i-1,j-1,k)+u(i-2,j-1,k))
c     &           /dx/cs(j-1)/dx/cs(j-1)
c     &           +(u(i-1,j,k)-2.*u(i-1,j-1,k)+u(i-1,j-2,k))/dy/dy
c            uafr10=(u(i+1,j-1,k)-2.*u(i,j-1,k)+u(i-1,j-1,k))
c     &           /dx/cs(j-1)/dx/cs(j-1)
c     &           +(u(i,j,k)-2.*u(i,j-1,k)+u(i,j-2,k))/dy/dy
c            uafr01=(u(i,j,k)-2.*u(i-1,j,k)+u(i-2,j,k))
c     &           /dx/cs(j)/dx/cs(j)
c     &           +(u(i-1,j+1,k)-2.*u(i-1,j,k)+u(i-1,j-1,k))/dy/dy
c            uafr11=(u(i+1,j,k)-2.*u(i,j,k)+u(i-1,j,k))
c     &           /dx/cs(j)/dx/cs(j)
c     &           +(u(i,j+1,k)-2.*u(i,j,k)+u(i,j-1,k))/dy/dy
c            vafr00=(v(i,j-1,k)-2.*v(i-1,j-1,k)+v(i-2,j-1,k))
c     &           /dx/cs(j-1)/dx/cs(j-1)
c     &           +(v(i-1,j,k)-2.*v(i-1,j-1,k)+v(i-1,j-2,k))/dy/dy
c            vafr10=(v(i+1,j-1,k)-2.*v(i,j-1,k)+v(i-1,j-1,k))
c     &           /dx/cs(j-1)/dx/cs(j-1)
c     &           +(v(i,j,k)-2.*v(i,j-1,k)+v(i,j-2,k))/dy/dy
c            vafr01=(v(i,j,k)-2.*v(i-1,j,k)+v(i-2,j,k))
c     &           /dx/cs(j)/dx/cs(j)
c     &           +(v(i-1,j+1,k)-2.*v(i-1,j,k)+v(i-1,j-1,k))/dy/dy
c            vafr11=(v(i+1,j,k)-2.*v(i,j,k)+v(i-1,j,k))
c     &           /dx/cs(j)/dx/cs(j)
c     &           +(v(i,j+1,k)-2.*v(i,j,k)+v(i,j-1,k))/dy/dy
c            visc_lp2=visc_lp2
c     &        +( (vafr11-vafr01)/dx/cs(j)/2.
c     &        +(vafr10-vafr00)/dx/cs(j-1)/2.
c     &        -(uafr11-uafr10)/dy/2.
c     &        -(uafr01-uafr00)/dy/2. )
c     &        *dzt(i,j,k)*faca*hduv
         enddo
      endif
      enddo
      enddo
      visc_lp=0.
      call sum_mst1(visc_lp,1)
c      call sum_mst(visc_lp2,1)
c
      btmstr=0.
      do i=3,iml-2
      do j=3,jml-2
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2.and.
     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2)then
         faca=0.
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        faca=faca+ex(i-1,j-1,1)*anhf(j-1)
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
     &        faca=faca+ex(i-1,j,1)*ashf(j)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        faca=faca+ex(i,j-1,1)*anhf(j-1)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
     &        faca=faca+ex(i,j,1)*ashf(j)
      do k=1,km
         vel00=0.
         vel10=0.
         vel01=0.
         vel11=0.
         if(exn(i-1,j-1).eq.k)vel00
     &        =dsqrt(u(i-1,j-1,k)**2+v(i-1,j-1,k)**2)*1.225d-3
         if(exn(i,j-1).eq.k)vel10
     &        =dsqrt(u(i,j-1,k)**2+v(i,j-1,k)**2)*1.225d-3
         if(exn(i-1,j).eq.k)vel01
     &        =dsqrt(u(i-1,j,k)**2+v(i-1,j,k)**2)*1.225d-3
         if(exn(i,j).eq.k)vel11
     &        =dsqrt(u(i,j,k)**2+v(i,j,k)**2)*1.225d-3
         btmstr=btmstr+(
     &        ( vel11*(v(i,j,k)*bcs+u(i,j,k)*sgn(j))
     &        -vel01*(v(i-1,j,k)*bcs+u(i-1,j,k)*sgn(j)) )
     &        /dx/cs(j)/2.
     &        +( vel10*(v(i,j-1,k)*bcs+u(i,j-1,k)*sgn(j-1))
     &        -vel00*(v(i-1,j-1,k)*bcs+u(i-1,j-1,k)*sgn(j-1)) )
     &        /dx/cs(j-1)/2.
     &        -( vel11*(u(i,j,k)*bcs-v(i,j,k)*sgn(j))
     &        -vel10*(u(i,j-1,k)*bcs-v(i,j-1,k)*sgn(j-1)) )
     &        /dy/2.
     &        -( vel01*(u(i-1,j,k)*bcs-v(i-1,j,k)*sgn(j))
     &        -vel00*(u(i-1,j-1,k)*bcs-v(i-1,j-1,k)*sgn(j)) )
     &        /dy/2.
     &        )*faca
      enddo
      endif
      enddo
      enddo
      call sum_mst(btmstr,1)
c
      zdwdz=0.
      fdwdz=0.
      do i=3,iml-2
      do j=3,jml-2
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2.and.
     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2)then
         faca=0.
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        faca=faca+ex(i-1,j-1,1)*anhf(j-1)
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
     &        faca=faca+ex(i-1,j,1)*ashf(j)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        faca=faca+ex(i,j-1,1)*anhf(j-1)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
     &        faca=faca+ex(i,j,1)*ashf(j)
         faca=faca/(areat(i,j,1)+1.-tex(i,j,1))
         do k=1,5
            zdwdz=zdwdz+zeta(i,j,k)*wl(i,j,1)/5.*faca
         enddo
         fdwdz=fdwdz+cor(j)*wl(i,j,1)*faca
      endif
      enddo
      enddo
      call sum_mst1(zdwdz,1)
      call sum_mst1(fdwdz,1)
c
      dwdxdvdz=0.
      dwdydudz=0.
      do i=3,iml-2
      do j=3,jml-2
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2-1.and.
     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2-1)then
      do k=2,km
         dwdxdvdz=dwdxdvdz-
     &        ( (wl(i+1,j+1,k)/(areat(i+1,j+1,k)+1.-tex(i+1,j+1,k))
     &        -wl(i,j+1,k)/(areat(i,j+1,k)+1.-tex(i,j+1,k)))
     &        /dx/cst(j+1)/2.
     &        +(wl(i+1,j,k)/(areat(i+1,j,k)+1.-tex(i+1,j,k))
     &        -wl(i,j,k)/(areat(i,j,k)+1.-tex(i,j,k)))
     &        /dx/cst(j)/2. )
     &        *(v(i,j,k-1)-v(i,j,k))*areauu(j)*ex(i,j,k)
         dwdydudz=dwdydudz+
     &        ( (wl(i+1,j+1,k)/(areat(i+1,j+1,k)+1.-tex(i+1,j+1,k))
     &        -wl(i+1,j,k)/(areat(i+1,j,k)+1.-tex(i+1,j,k)))
     &        /dy/2.
     &        +(wl(i,j+1,k)/(areat(i,j+1,k)+1.-tex(i,j+1,k))
     &        -wl(i,j,k)/(areat(i,j,k)+1.-tex(i,j,k)))
     &        /dy/2. )
     &        *(u(i,j,k-1)-u(i,j,k))*areauu(j)*ex(i,j,k)
      enddo
      endif
      enddo
      enddo
      call sum_mst1(dwdxdvdz,1)
      call sum_mst1(dwdydudz,1)

      call wait_t3d_w1(zeta2x,11)
      call wait_t3d_e1(zeta2x,12)
      call wait_t3d_s1(zeta2y,13)
      call wait_t3d_n1(zeta2y,14)
c
      visc_bh=0.
      do i=3,iml-2
      do j=3,jml-2
      if(lw(ip_x)+i-1.ge.ii1.and.lw(ip_x)+i-1.le.ii2.and.
     &        ls(ip_y)+j-1.ge.jj1.and.ls(ip_y)+j-1.le.jj2)then
c         faca=0.
c         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
c     &        faca=faca+ex(i-1,j-1,1)*anhf(j-1)
c         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
c     &        faca=faca+ex(i-1,j,1)*ashf(j)
c         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
c     &        faca=faca+ex(i,j-1,1)*anhf(j-1)
c         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
c     &        faca=faca+ex(i,j,1)*ashf(j)
         do k=1,km
c
         facv=0.
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        facv=facv+ex(i-1,j-1,k)*anhf(j-1)*dzu(i-1,j-1,k)
         if(lw(ip_x)+i-1.ge.ii1+1.and.ls(ip_y)+j-1.le.jj2-1)
     &        facv=facv+ex(i-1,j,k)*ashf(j)*dzu(i-1,j,k)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.ge.jj1+1)
     &        facv=facv+ex(i,j-1,k)*anhf(j-1)*dzu(i,j-1,k)
         if(lw(ip_x)+i-1.le.ii2-1.and.ls(ip_y)+j-1.le.jj2-1)
     &        facv=facv+ex(i,j,k)*ashf(j)*dzu(i,j,k)
c
            visc_bh=visc_bh-hduv_bh_x*(
     &        ( (zeta2x(i+1,j,k)-zeta2x(i,j,k))*tex(i+1,j,k)
     &        -(zeta2x(i,j,k)-zeta2x(i-1,j,k))*tex(i-1,j,k) )
     &        /dx/cst(j)/dx/cst(j)
     &        +( (zeta2y(i,j+1,k)-zeta2y(i,j,k))*tex(i,j+1,k)
     &        -(zeta2y(i,j,k)-zeta2y(i,j-1,k))*tex(i,j-1,k) )
     &           /dy/dy )*facv
         enddo
         endif
      enddo
      enddo
      call sum_mst1(visc_bh,1)
c
      adv=adv_w+adv_e+adv_s+adv_n
      betav2=betav2w+betav2e+betav2s+betav2n
      allterms=adv_w+adv_e+adv_s+adv_n+curltau+betav2+btmstr
     &     +visc_lp+visc_bh+zdwdz+fdwdz+dwdxdvdz+dwdydudz
c
c
      if(ip.eq.imaster)then
         write(*,*)' '
         write(*,*)nkai
         write(*,*)'circulation:       ',circ
         write(*,*)'adv:               ',adv,
     &        ': ',adv_w,adv_e,adv_s,adv_n
c         write(*,*)'betav:             ',betav,betav2
         write(*,*)'betav:             ',betav2,
     &        ': ',betav2w,betav2e,betav2s,betav2n
         write(*,*)'curltau,bottom:    ',curltau,btmstr
         write(*,*)'viscosity lp,bh:   ',visc_lp,visc_bh
c         write(*,*)'visc lp 2:         ',visc_lp2
         write(*,*)'zdwdz,fdwdz:       ',zdwdz,fdwdz
         write(*,*)'dwdxdvdz,dwdydudz: ',dwdxdvdz,dwdydudz
         write(*,*)'total:             ',allterms
         write(vortegf)ahour,adv,adv_w,adv_e,adv_s,adv_n,
     &        betav2,curltau,btmstr,visc_lp,visc_bh,
     &        zdwdz,fdwdz,dwdxdvdz,dwdydudz,allterms
      endif
c
c
      endif
c
      return
      end
#endif
c
c
