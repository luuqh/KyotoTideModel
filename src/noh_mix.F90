ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     Noh-mixed layer model                                         c
c                                                                   c
c     $Id: noh_mix.F90 2 2007-11-15 13:07:08Z ishikawa $
c                                                                   c
c     tke is loop only this subroutine                              c
c     vdts, vduv is output for tracli.F                             c
c     request for tracli  rhoo,rhodh,rhouh                          c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine noh_mix
      
      use param
      use mod_mpi
      implicit none

#include "common.h"
   
      integer i,j,k,n,modmtnoh,matsunoh
      real*8 sst,drdn,drup,dpth,wtx,wty,
     &     disp,sprod,drho,fact1,tkevt,facup,facdn,vdf2,dudz,dvdz,
     &     facw,face,facs,facn
#ifdef NOHOLD
      real*8 ri
#else
      real*8 ri(im,jm,0:km+1)
#endif
      real*8 vddsf,vddtc,vddsc,ds_temp
      integer jsn,jen
      real*8 tkecoef
      integer ntkeloop
      real*8 tkem(im,jm,0:km+1),tkemb(im,jm,0:km+1)!,dttke,c2dttke
c
ccc   for implicit vertical diffusion
      real*8 tria3(im,km),trib3(im,km),tric3(im,km),
     &  tribet3(im),trigam3(im,km)
c
      call exch_t2d_s1p(wsx,wsy,5)
c
      if(matsno.eq.1.and.matsn2.eq.0) then
         do k = 1,km
         do j = 1,jml
         do i = 1,iml
            tkeb(i,j,k) = tke(i,j,k)
         enddo
         enddo
         enddo
      endif
c
      do k = 0,km+1
      do j = 1,jm
      do i = 1,im
         tkem(i,j,k) = dmax1(tke(i,j,k),tkemin)
         tkemb(i,j,k) = dmax1(tke(i,j,k),tkemin)
      enddo
      enddo
      enddo
      
      call wait_t2d_s1p(wsx,wsy,5)
      call exch_t2d_w1p(wsx,wsy,7)

      do k = 0,km+1
      do j = 1,jm
      do i = 1,im
         tkea(i,j,k)=0.d0
         tked(i,j,k)=0.d0
         vdts(i,j,k)=0.d0
         vdtke(i,j,k)=0.d0
         vddt(i,j,k)=0.d0
         vdds(i,j,k)=0.d0       
#ifndef NOHOLD
         ri(i,j,k)=0.d0
#endif
      enddo
      enddo
      enddo
c
ccc   compute mixed layer depth
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,1) .eq. 1.d0) then
         sst = dmax1(t(i,j,1),-3.d0)
         sst = dmin1(sst,30.d0)
         rhomix(i,j) = 1.d0-tex(i,j,1)+tex(i,j,1)*(
     &        rhoo(i,j,kmlmin)+drfit(1) +drfit(2)*sst
     &        +drfit(3)*sst*sst+drfit(4)*sst**3
     &        +drfit(5)*sst**4 )
c
      do k = km,kmlmin,-1
         drdn = rhomix(i,j) - rhoo(i,j,k)
         drup = rhomix(i,j) - rhoo(i,j,k-1)
c
         if(tex(i,j,k) .eq. 0.) then
           depml(i,j)=dep(k-1)+hcl(i,j)+0.5d0*dzt(i,j,k-1)
         endif
         
         if(drdn.lt.0.d0 .and. drup.gt.0.d0)then
            if(k-1.le.kadp) then
               dpth=dep(k-1)*(1.d0+hcl(i,j)/depadp)+0.5d0*dzt(i,j,k-1)
            else
               dpth=dep(k-1)+hcl(i,j)+0.5d0*dzt(i,j,k-1)
            endif
            depml(i,j) = dpth + 0.5d0*(dzt(i,j,k-1)+dzt(i,j,k))
     &           *drup/(rhoo(i,j,k)-rhoo(i,j,k-1))
         endif
      enddo
c
         depml(i,j)=dmax1(depml(i,j),10.d2)
c         depml(i,j)=dmin1(depml(i,j),2000.d2)
         depml(i,j)=dmin1(depml(i,j),dep(km+1))
      endif
      enddo
      enddo
c
#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do k=1,km
      do j=1,jml
      do i=1,indg+2
         tkeawbt(i,j,k) = tkewbc(i,j,k,mb1)
     &        + (tkewbc(i,j,k,mb2)-tkewbc(i,j,k,mb1)) * cbf
      enddo
      enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do k=1,km
      do j=1,jml
      do i=1,indg+2
         tkeaebt(i,j,k) = tkeebc(i,j,k,mb1)
     &        + (tkeebc(i,j,k,mb2)-tkeebc(i,j,k,mb1)) * cbf
      enddo
      enddo
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do k=1,km
      do i=1,iml
      do j=1,jndg+2
        tkeasbt(i,j,k) = tkesbc(i,j,k,mb1)
     &           + (tkesbc(i,j,k,mb2)-tkesbc(i,j,k,mb1)) * cbf
      enddo
      enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do k=1,km
      do i=1,iml
      do j=1,jndg+2
         tkeanbt(i,j,k) = tkenbc(i,j,k,mb1)
     &        + (tkenbc(i,j,k,mb2)-tkenbc(i,j,k,mb1)) * cbf
      enddo
      enddo
      enddo
      endif
#endif

      call wait_t2d_w1p(wsx,wsy,7)

      c2dttke=dttke*2.d0
      modmtnoh=5
      ntkeloop=max(int(c2dtts/c2dttke),1)
      do 1201 n=1,ntkeloop

#ifdef NESTED
      tkecoef=dble(n)/dble(ntkeloop)
#endif
c
      matsunoh=0
      if(mod(n,modmtnoh).eq.0)then
         matsunoh=1
         c2dttke=dttke
         do i=2,iml-1
         do j=2,jml-1
         do k=1,km
            tkemb(i,j,k)=tkem(i,j,k)
         enddo
         enddo
         enddo
      endif
 1301 continue
c
      call exch_t3d_n1p(tkem,tkemb,3,4)
c
ccc   surface vdts/vdtke
      do j = 3,jml-2
      do i = 3,iml-2
         dleng(i,j,1) = ckrmn*(0.25*dzt(i,j,1)+rough0)
     &        /(1.d0+ckrmn*(0.25*dzt(i,j,1)+rough0)
     &        /(depml(i,j)+1.-tex(i,j,1)))
c
#ifdef NOHOLD
         vdts(i,j,1)=s0/cpr*dleng(i,j,1)*dsqrt(2.d0*tkem(i,j,1))
#else
         vdts(i,j,1)=s0/pr0*dleng(i,j,1)*dsqrt(2.d0*tkem(i,j,1))
#endif
         vdtke(i,j,1)=s0/sigma*dleng(i,j,1)*dsqrt(2.d0*tkem(i,j,1))
         vdts(i,j,1)=dmax1(vdts(i,j,1),vdtsmn)
         vdtke(i,j,1)=dmax1(vdtke(i,j,1),vdtsmn)
         vdts(i,j,1)=tex(i,j,1)*dmin1(vdts(i,j,1),vdtsmx)
         vdtke(i,j,1)=tex(i,j,1)*dmin1(vdtke(i,j,1),vdtsmx)
      enddo
      enddo

      call wait_t3d_n1p(tkem,tkemb,3,4)
      call exch_t3d_e1p(tkem,tkemb,8,9)

      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,1).eq.1.d0)then
ccc   surface condition of tke
         wtx = (anhf(j-1)*(
     &     ex(i-1,j-1,1)*wsx(i-1,j-1)**2+ex(i,j-1,1)*wsx(i,j-1)**2)
     &     +ashf(j)*(ex(i,j,1)*wsx(i,j)**2+ex(i-1,j,1)*wsx(i-1,j)**2)
     &     )/(areat(i,j,1)+1.-tex(i,j,1))
         wty = (anhf(j-1)*(
     &     ex(i-1,j-1,1)*wsy(i-1,j-1)**2+ex(i,j-1,1)*wsy(i,j-1)**2)
     &     +ashf(j)*(ex(i,j,1)*wsy(i,j)**2+ex(i-1,j,1)*wsy(i-1,j)**2)
     &     )/(areat(i,j,1)+1.-tex(i,j,1))
c
#ifdef NOHOLD
         disp = c0/(dleng(i,j,1)+1.-tex(i,j,1))
     &        *dsqrt(2.d0*tkem(i,j,1))*2.d0*tkemb(i,j,1)*0.2d0
         sprod = (wtx+wty)*ddmnar**2
     &        /(cpr*vdts(i,j,1)+1.d0-tex(i,j,1))
#else
         disp = c0_e/dleng_e
     &        *dsqrt(2.d0*tkem(i,j,1))*2.d0*tkemb(i,j,1)
         sprod = (wtx+wty)*ddmnar**2
     &        /(pr0*vdts(i,j,1)+1.d0-tex(i,j,1))
#endif
c
         tkea(i,j,1) = tkea(i,j,1)+areat(i,j,1)*(cm0*fricv3(i,j)
     &        + 0.5*dzt(i,j,1)*(-disp+sprod) )
      endif
c
ccc   double diffusion
         vddt(i,j,1)=vdts(i,j,1)*tex(i,j,1)
         vdds(i,j,1)=vdts(i,j,1)*tex(i,j,1)
      enddo
      enddo
c
ccc   vdts/vdtke(2:km)
      do k = 2,km
      do j = 3,jml-2
      do i = 3,iml-2
      drho = tex(i,j,k)*(rhodh(i,j,k-1)-rhouh(i,j,k))*dzzr(k)

      if(k.le.kadp) then
         dpth = dep(k)*(1.+hcl(i,j)/depadp)
      else
         dpth = dep(k) + hcl(i,j)
      endif
      dleng(i,j,k) = ckrmn*(dpth+rough0)
     &     /(1.d0+ckrmn*(dpth+rough0)/
     &     (depml(i,j)+1.-tex(i,j,k)))
      if(dep(k).gt.depml(i,j))dleng(i,j,k)=100.d2
c
#ifdef NOHOLD
      ri = -tex(i,j,k)*0.5d0*grav*drho*dleng(i,j,k)**2
     &     /(dmax1(tkem(i,j,k),tkemin)+1.-tex(i,j,k))
      ri=dmax1(ri,0.d0)
      dleng(i,j,k)=dleng(i,j,k)/dsqrt(1.d0+calph*ri)
c
      vdts(i,j,k)=s0/cpr*dleng(i,j,k)*sqrt(2.d0*tkem(i,j,k))
      vdtke(i,j,k)=s0/sigma*dleng(i,j,k)*sqrt(2.d0*tkem(i,j,k))
      disp=tex(i,j,k)*c0/(dleng(i,j,k)+1.-tex(i,j,k))
     &     *dsqrt(2.d0*tkem(i,j,k))*2.d0*tkemb(i,j,k)
#else
      ri(i,j,k)=-tex(i,j,k)*0.5d0*grav*drho*dleng(i,j,k)**2
     &     /(dmax1(tkem(i,j,k),tkemin)+1.-tex(i,j,k))
      ri(i,j,k)=dmax1(ri(i,j,k),0.d0)
      prandtl=pr0*dsqrt(1.d0+calph_p*ri(i,j,k))
c
      vdts(i,j,k)=s0/prandtl*dleng(i,j,k)/dsqrt(1.d0+calph_b*ri(i,j,k))
     &     *sqrt(2.d0*dmax1(tkem(i,j,k),tkemin))
      vdtke(i,j,k)=s0/sigma*dleng(i,j,k)/dsqrt(1.d0+calph_b*ri(i,j,k))
     &     *sqrt(2.d0*dmax1(tkem(i,j,k),tkemin))
      disp=tex(i,j,k)*c0_e/dleng_e*dsqrt(1.d0+calph_e*ri(i,j,k))
     &     *dsqrt(2.d0*dmax1(tkem(i,j,k),tkemin))*2.d0*tkemb(i,j,k)
#endif
c
      vdts(i,j,k)=dmax1(vdts(i,j,k),vdtsmn)
      vdtke(i,j,k)=dmax1(vdtke(i,j,k),vdtsmn)
      vdts(i,j,k)=tex(i,j,k)*dmin1(vdts(i,j,k),vdtsmx)
      vdtke(i,j,k)=tex(i,j,k)*dmin1(vdtke(i,j,k),vdtsmx)

ccc   double diffusion parameterization
      vddt(i,j,k)=vdts(i,j,k)
      vdds(i,j,k)=vdts(i,j,k)
!      dratio=tex(i,j,k)*expt/exps*(t(i,j,k-1)-t(i,j,k))
!     &     /dmax1((s(i,j,k-1)-s(i,j,k)),1.d-5)
       ds_temp=s(i,j,k-1)-s(i,j,k)
       if(ds_temp==0.) ds_temp=1.d-10
       dratio=tex(i,j,k)*expt/exps*(t(i,j,k-1)-t(i,j,k))/ds_temp
      if(dratio.gt.1..and.dratio.lt.dratioc.and.
     &     t(i,j,k-1).gt.t(i,j,k).and.s(i,j,k-1).gt.s(i,j,k))then
         vddsf=vddsmx*(1.d0-((dratio-1.d0)/(dratioc-1.d0))**2)**3
         vdds(i,j,k)=dmax1(vdds(i,j,k),vddsf)*tex(i,j,k)
         vddt(i,j,k)=dmax1(vddt(i,j,k),vddsf*0.7)*tex(i,j,k)
      else if(dratio.gt.0..and.dratio.lt.1..and.
     &     t(i,j,k-1).lt.t(i,j,k).and.s(i,j,k-1).lt.s(i,j,k))then
         vddtc=0.909*vddtmn*exp(4.6*exp(-0.54*(1./dratio-1.)))
         vddt(i,j,k)=dmax1(vddt(i,j,k),vddtc)*tex(i,j,k)
         if(dratio.ge.0.5.and.dratio.lt.1.)then
            vddsc=vddtc*(1.85-0.85/dratio)*dratio
         else if(dratio.gt.0..and.dratio.lt.0.5)then
            vddsc=vddtc*0.15*dratio
         endif
         vdds(i,j,k)=dmax1(vdds(i,j,k),vddsc)*tex(i,j,k)
      endif
      if(vddt(i,j,k).gt.vdts(i,j,k).or.
     &     vdds(i,j,k).gt.vdts(i,j,k))then
         vdts(i,j,k)=dmax1(vddt(i,j,k),vdds(i,j,k))
#ifdef NOHOLD
         vdtke(i,j,k)=dmax1(vdts(i,j,k)*cpr/sigma,vdtsmn)
#else
         vdtke(i,j,k)=dmax1(vdts(i,j,k)*prandtl/sigma,vdtsmn)
#endif
         vdtke(i,j,k)=tex(i,j,k)*dmin1(vdtke(i,j,k),vdtsmx)
      endif

ccc   tidal mixing in cril
      vdts(i,j,k)=dmax1(vdts(i,j,k),vdtide(i,j,k))*tex(i,j,k)


ccc   tsujino shceme
c      vdts(i,j,k)=tex(i,j,k)*dmax1(vdts(i,j,k),vdtsuji(k))
c      vddt(i,j,k)=tex(i,j,k)*dmax1(vddt(i,j,k),vdtsuji(k))
c      vdds(i,j,k)=tex(i,j,k)*dmax1(vdds(i,j,k),vdtsuji(k))

      tkea(i,j,k) = tkea(i,j,k) + (grav*vdts(i,j,k)*drho-disp)
     &  *0.5*(areat(i,j,k-1)*dzt(i,j,k-1)+areat(i,j,k)*dzt(i,j,k))

#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59 .and.
     &   i==33 .and. j==51 ) then
           write(*,*) nkai,'vdts',k,vdts(i,j,k),tkem(i,j,k),drho,
     &      ri(i,j,k)
      endif
#endif


      enddo
      enddo
      enddo
#ifdef DEBUG1
      if(nkai>482095 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 1,km
          write(*,*) 'tkea-local',tkea(i,j,k),vdts(i,j,k)
        enddo
        write(*,*)
      endif
#endif
      call wait_t3d_e1p(tkem,tkemb,8,9)

ccc   horizontal advection and diffusion of tke
      do k = 2,km
      do j = 3,jml-2
      do i = 3,iml-2
         if(ustar(i,j,k-1)+ustar(i,j,k).gt.0.d0) then
            facw = (1.d0+hupp)*0.5d0
            face = (1.d0-hupp)*0.5d0
         else
            facw = (1.d0-hupp)*0.5d0
            face = (1.d0+hupp)*0.5d0
         endif

         eust(i,j,k)=0.25d0*(
     &        (ustar(i,j,k-1)*(ex(i,j,k-1)+ex(i,j-1,k-1))
     &        +ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k)))
     &        *(facw*tkem(i,j,k)+face*tkem(i+1,j,k))
     &        -dydxr*hdts*cstr(j)
     &        *(dzu(i,j,k)+dzu(i,j-1,k)+dzu(i,j,k-1)+dzu(i,j-1,k-1))
     &        *(tkemb(i+1,j,k)-tkemb(i,j,k)) )

         if(vstar(i,j,k-1)+vstar(i,j,k).gt.0.d0) then
            facs = (1.d0+hupp)*0.5d0
            facn = (1.d0-hupp)*0.5d0
         else
            facs = (1.d0-hupp)*0.5d0
            facn = (1.d0+hupp)*0.5d0
         endif

         evst(i,j,k)=0.25d0*(
     &        (vstar(i,j,k-1)*(ex(i,j,k-1)+ex(i-1,j,k-1))
     &        +vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k)))
     &        *(facs*tkem(i,j,k)+facn*tkem(i,j+1,k))
     &        -dydxr*hdts*cstr(j)
     &        *(dzu(i,j,k)+dzu(i-1,j,k)+dzu(i,j,k-1)+dzu(i-1,j,k-1))
     &        *(tkemb(i,j+1,k)-tkemb(i,j,k)) )

      enddo
      enddo
      enddo

      call exch_t3d_s1(evst,2)
      call exch_t3d_w1(eust,10)

ccc   volume change of tke box
      do k = 1,kadp+1
      do j = 3,jml-2
      do i = 3,iml-2
         if(k.eq.kadp+1) then
            fact1 = (dz(k-1))/depadp*0.5d0
         else
            fact1 = (dz(k-1)+dz(k))/depadp*0.5d0
         endif
         tkea(i,j,k)=tkea(i,j,k)
     &        -tkemb(i,j,k)*(wl(i,j,1)+wflux(i,j)*areat(i,j,1))
     &        *fact1*tex(i,j,k)
      enddo
      enddo
      enddo
#ifdef DEBUG1
      if(nkai>482095 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 1,km
          write(*,*) 'tkea-change',tkea(i,j,k)
        enddo
        write(*,*)
      endif
#endif
ccc   vertical advection and diffusion
      do k = 1,km-1
      do j = 3,jml-2
      do i = 3,iml-2
         if((wl(i,j,k)+wl(i,j,k+1)).gt.0.d0) then
            facup = (1.d0-vupp)*0.5d0
            facdn = (1.d0+vupp)*0.5d0
         else
            facup = (1.d0+vupp)*0.5d0
            facdn = (1.d0-vupp)*0.5d0
         endif

         vdf2=0.5*(vdtke(i,j,k)+vdtke(i,j,k+1))
     $        *areat(i,j,k)/(dzt(i,j,k)+1.-tex(i,j,k))
         if(k.eq.1) then
            tkevt=0.5*wl(i,j,k+1)
     &           *(facup*tkem(i,j,k)+facdn*tkem(i,j,k+1))
     &           -(tkemb(i,j,k)-tkemb(i,j,k+1))*vdf2
         else
            tkevt=0.5*(wl(i,j,k)+wl(i,j,k+1))
     &           *(facup*tkem(i,j,k)+facdn*tkem(i,j,k+1))
     &           -(tkemb(i,j,k)-tkemb(i,j,k+1))*vdf2
         endif
         tkea(i,j,k)=tkea(i,j,k)+tex(i,j,k)*tkevt
         tkea(i,j,k+1)=tkea(i,j,k+1)-tex(i,j,k)*tkevt

      enddo
      enddo
      enddo
#ifdef DEBUG1
      if(nkai>482095 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 1,km
          write(*,*) 'tkea-vt',tkea(i,j,k)
        enddo
        write(*,*)
      endif
#endif
      call wait_t3d_s1(evst,2)
      call wait_t3d_w1(eust,10)
c
ccc   vartical shear production
      do k = 2,km
      do j = 3,jml-2
      do i = 3,iml-2
         dudz = anhf(j-1)*(
     &        ex(i-1,j-1,k)*(u(i-1,j-1,k-1)-u(i-1,j-1,k))**2
     &        +ex(i,j-1,k)*(u(i,j-1,k-1)-u(i,j-1,k))**2)
     &        +ashf(j)*(ex(i,j,k)*(u(i,j,k-1)-u(i,j,k))**2
     &        +ex(i-1,j,k)*(u(i-1,j,k-1)-u(i-1,j,k))**2)
         dvdz = anhf(j-1)*(
     &        ex(i-1,j-1,k)*(v(i-1,j-1,k-1)-v(i-1,j-1,k))**2
     &        +ex(i,j-1,k)*(v(i,j-1,k-1)-v(i,j-1,k))**2)
     &        +ashf(j)*(ex(i,j,k)*(v(i,j,k-1)-v(i,j,k))**2
     &        +ex(i-1,j,k)*(v(i-1,j,k-1)-v(i-1,j,k))**2)
#ifdef NOHOLD
         sprod = cpr*vdts(i,j,k)*(dudz+dvdz)*dzzr(k)
#else
         prandtl=pr0*dsqrt(1.d0+calph_p*ri(i,j,k))
         sprod = prandtl*vdts(i,j,k)*(dudz+dvdz)*dzzr(k)
#endif
         tkea(i,j,k) = tkea(i,j,k) + sprod
     &     +(eust(i-1,j,k)-eust(i,j,k)+evst(i,j-1,k)-evst(i,j,k))
      enddo
      enddo
      enddo
      
#ifdef DEBUG1
      if(nkai>482095 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 1,km
          write(*,*) 'tkea-last',tkea(i,j,k),eust(i-1,j,k),
     &     eust(i,j,k),evst(i,j-1,k),evst(i,j,k)
        enddo
        write(*,*)
      endif
#endif

#ifdef NESTED
! south
      if(ip_y .eq. 0) then
      do k=1,km
      do i=3,iml-2
      do j=3,jndg+2
            tkea(i,j+2,k) = tkea(i,j+2,k) + 
     &           tex(i,j+2,k)*chfq*(tkesbt(i,j,k)-tkem(i,j+2,k))
     &           *areat(i,j+2,k)*dzt(i,j+2,k)
      enddo
      enddo
      enddo      
! avoid doublt counting
        jsn = jndg+5
      else
        jsn = 3
      endif
! north
      if(ip_y .eq. jpe-1) then
      do k=1,km
      do i=3,iml-2
      do j=3,jndg+2
            tkea(i,jml-1-j,k) = tkea(i,jml-1-j,k) + 
     &        tex(i,jml-1-j,k)*chfq*(tkenbt(i,j,k)-tkem(i,jml-1-j,k))
     &       *areat(i,jml-1-j,k)*dzt(i,jml-1-j,k)

      enddo
      enddo
      enddo
! avoid doublt counting
        jen = jml-jndg-4
      else
        jen = jml-2
      endif
! west      
      if(ip_x .eq. 0) then
      do k=1,km
      do j = jsn,jen
      do i=3,indg+2
            tkea(i+2,j,k) = tkea(i+2,j,k) + 
     &           tex(i+2,j,k)*chfq*(tkewbt(i,j,k)-tkem(i+2,j,k))
     &           *areat(i+2,j,k)*dzt(i+2,j,k)
      enddo
      enddo
      enddo      
      endif
! east
      if(ip_x .eq. ipe-1) then
      do k=1,km
      do j = jsn,jen
      do i=3,indg+2
            tkea(iml-1-i,j,k) = tkea(iml-1-i,j,k) + 
     &       tex(iml-1-i,j,k)*chfq*(tkeebt(i,j,k)-tkem(iml-1-i,j,k))
     &      *areat(iml-1-i,j,k)*dzt(iml-1-i,j,k)
      enddo
      enddo
      enddo
      endif
#endif


ccc   solving tri-diag matrix inversion for implicit vertical mixing
c
      do j = 3,jml-2
ccc   form matrix
c
ccc   k = 1
      do i = 3,iml-2
         tric3(i,1) = -(vdtke(i,j,1)+vdtke(i,j,2))*areat(i,j,1)
     &     /(dzt(i,j,1)*(volt(i,j,1)+volt(i,j,2))+1.-tex(i,j,1))
#ifdef NOHOLD
         trib3(i,1) = tex(i,j,1)/c2dttke-tric3(i,1)
     &        +0.06d0*dsqrt(2.d0*tkem(i,j,1))*2.d0
     &        /(dleng(i,j,1)+1.-tex(i,j,1))*0.2d0
#else
         trib3(i,1) = tex(i,j,1)/c2dttke-tric3(i,1)
     &        +c0_e*dsqrt(2.d0*tkem(i,j,1))*2.d0
     &        /dleng_e
#endif
      enddo
c
      do k = 2,km-1
      do i = 3,iml-2
         tria3(i,k) = -tex(i,j,k)*(vdtke(i,j,k-1)+vdtke(i,j,k))
     &     *areat(i,j,k)/(dzt(i,j,k-1)*(volt(i,j,k-1)+volt(i,j,k))
     &        +1.-tex(i,j,k-1))
         tric3(i,k) = -tex(i,j,k+1)*(vdtke(i,j,k)+vdtke(i,j,k+1))
     &     *areat(i,j,k+1)/(dzt(i,j,k)*(volt(i,j,k)+volt(i,j,k+1))
     &        +1.-tex(i,j,k))
#ifdef NOHOLD
         trib3(i,k) = tex(i,j,k)/c2dttke-tria3(i,k)-tric3(i,k)
     &        +0.06d0*dsqrt(2.d0*tkem(i,j,k))*2.d0
     &        /(dleng(i,j,k)+1.-tex(i,j,k))
#else
         trib3(i,k) = tex(i,j,k)/c2dttke-tria3(i,k)-tric3(i,k)
     &        +c0_e*dsqrt(2.d0*tkem(i,j,k))*2.d0
     &        /dleng_e*dsqrt(1.d0+calph_e*ri(i,j,k))
#endif
      enddo
      enddo
c
ccc   k = km
      do i = 3,iml-2
         tria3(i,km) = -tex(i,j,km)*(vdtke(i,j,km-1)+vdtke(i,j,km))
     &  *areat(i,j,km)/(dzt(i,j,km-1)*(volt(i,j,km-1)+volt(i,j,km))
     &        +1.-tex(i,j,km-1))
#ifdef NOHOLD
         trib3(i,km) = tex(i,j,km)/c2dttke-tria3(i,km)
     &        +0.06d0*dsqrt(2.d0*tkem(i,j,km))*2.d0
     &        /(dleng(i,j,km)+1.-tex(i,j,km))
#else
         trib3(i,km) = tex(i,j,km)/c2dttke-tria3(i,km)
     &        +c0_e*dsqrt(2.d0*tkem(i,j,km))*2.d0
     &        /dleng_e*dsqrt(1.d0+calph_e*ri(i,j,km))
#endif
      enddo


ccc   solve equation

      do i = 3,iml-2
      if(tex(i,j,1).eq.1.d0) then
         tkea(i,j,1) = tkea(i,j,1)*2.*voltr(i,j,1)

         tribet3(i) = trib3(i,1)
         tked(i,j,1) = tkea(i,j,1)/tribet3(i)
      endif
      enddo

      do k = 2,km
      do i = 3,iml-2
      if(tex(i,j,k).eq.1.d0) then
         tkea(i,j,k) = tkea(i,j,k)*2.
     &        /(volt(i,j,k-1)+volt(i,j,k)+1.-tex(i,j,k-1))

         trigam3(i,k) = tric3(i,k-1)/tribet3(i)
         tribet3(i) = trib3(i,k) - tria3(i,k)*trigam3(i,k)
         tked(i,j,k)=(tkea(i,j,k)-tria3(i,k)*tked(i,j,k-1))/tribet3(i)
      endif
      enddo
      enddo

      do k = km-1,1,-1
      do i = 3,iml-2
      if(tex(i,j,k+1).eq.1.d0) then
         tked(i,j,k)=tked(i,j,k)-trigam3(i,k+1)*tked(i,j,k+1)
      endif
      enddo
      enddo

      enddo

      do j = 3,jml-2
      do i = 3,iml-2
         tkea(i,j,1)=dmax1(tkemb(i,j,1)+tked(i,j,1),tkemin)
         tkea(i,j,1)=tex(i,j,1)*dmin1(tkea(i,j,1),1.d6)
      enddo
      enddo

      do k = 2,km
      do j = 3,jml-2
      do i = 3,iml-2
         tkea(i,j,k)=dmax1(tkemb(i,j,k)+tked(i,j,k),tkemin)
         tkea(i,j,k)=tex(i,j,k)*dmin1(tkea(i,j,k),1.d6)
      enddo
      enddo
      enddo

#ifdef DEBUG1
      if(nkai>482095 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 1,km
          write(*,*) 'tke-last',tkea(i,j,k),tkemb(i,j,k),tked(i,j,k)
        enddo
        write(*,*)
      endif
#endif

#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do k=1,km
      do j=3,jml-2
      do i=1,2
      tkea(i+2,j,k) = (tkewbt(i,j,k)+(tkeawbt(i,j,k)-tkewbt(i,j,k))
     $        *tkecoef)*tex(i+2,j,k)
      enddo
      enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do k=1,km
      do j=3,jml-2
      do i=1,2
      tkea(iml-1-i,j,k) = (tkeebt(i,j,k)+(tkeaebt(i,j,k)-tkeebt(i,j,k))
     $        *tkecoef)*tex(iml-1-i,j,k)
      enddo
      enddo
      enddo      
      endif
!  south
      if(ip_y .eq. 0) then
      do k = 1,km
      do i = 1,iml
      do j=1,2
      tkea(i,j+2,k) = (tkesbt(i,j,k)+(tkeasbt(i,j,k)-tkesbt(i,j,k))
     $        *tkecoef)*tex(i,j+2,k)
      enddo
      enddo
      enddo      
      endif
! north
      if(ip_y .eq. jpe-1) then
      do k = 1,km
      do i = 1,iml
      do j=1,2
      tkea(i,jml-1-j,k) = (tkenbt(i,j,k)+(tkeanbt(i,j,k)-tkenbt(i,j,k))
     $        *tkecoef)*tex(i,jml-1-j,k)
      enddo
      enddo
      enddo
      
      endif
#endif

      do k=1,km
      do j = 1,jml
      do i = 1,iml
        if(matsunoh.eq.1) then
          tkem(i,j,k)=dmax1(tkea(i,j,k),tkemin)        
        else
          tkemb(i,j,k)=tkem(i,j,k)
          tkem(i,j,k)=dmax1(tkea(i,j,k),tkemin)
        endif
      enddo
      enddo
      enddo
#ifdef DEBUG1
      if(nkai>482101 .and. nkai < 482110) then
!      if(nkai>482016 .and. nkai < 482112) then
        if(ip==59) then
          write(90) nkai,n
          write(90) tkem
        endif
      endif
#endif

      if(matsunoh.eq.1)then
         matsunoh=0
         c2dttke=dttke*2.d0
         goto 1301
      endif

 1201 continue

      call exch_t3d_n1(vdts,1)
      call wait_t3d_n1(vdts,1)
      call exch_t3d_e1(vdts,11)
#ifndef NOHOLD
      call exch_t3d_n1(ri,6)
      call wait_t3d_n1(ri,6)
      call exch_t3d_e1(ri,12)
#endif
c
ccc   transport things to prepare next step
      if(matsno.eq.1.or.matsn2.eq.1) then
         do k = 1,km
         do j = 1,jml
         do i = 1,iml
            tke(i,j,k)=dmax1(tkea(i,j,k),tkemin)
         enddo
         enddo
         enddo
      else 
         do k = 1,km
         do j = 1,jml
         do i = 1,iml
            tkeb(i,j,k)=tke(i,j,k)
            tke(i,j,k)=dmax1(tkea(i,j,k),tkemin)
         enddo
         enddo
         enddo
      endif
c
#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do k=1,km
      do j=1,jm
      do i=1,indg+2
         tkewbt(i,j,k) = tkeawbt(i,j,k)
      enddo
      enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do k=1,km
      do j=1,jml
      do i=1,indg+2
         tkeebt(i,j,k) = tkeaebt(i,j,k)
      enddo
      enddo
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do k=1,km
      do i=1,iml
      do j=1,jndg+2
        tkesbt(i,j,k) = tkeasbt(i,j,k)
      enddo
      enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do k=1,km
      do i=1,iml
      do j=1,jndg+2
            tkenbt(i,j,k) = tkeanbt(i,j,k)
      enddo
      enddo
      enddo
      endif
#endif
c
ccc   vduv
      call wait_t3d_e1(vdts,11)
#ifndef NOHOLD
      call wait_t3d_e1(ri,12)
#endif
c
      do k = 2,km
      do j = 3,jml-2
      do i = 3,iml-2
#ifdef NOHOLD
         vduv(i,j,k)=ex(i,j,k)/areauu(j)*cpr*(
     &        ashf(j)*(vdts(i,j,k)+vdts(i+1,j,k))
     &        +anhf(j)*(vdts(i,j+1,k)+vdts(i+1,j+1,k)))
#else
         vduv(i,j,k)=0.d0
         prandtl=pr0*dsqrt(1.d0+calph_p*ri(i,j,k))
         vduv(i,j,k)=ex(i,j,k)/areauu(j)*(
     &     ashf(j)*(vdts(i,j,k)*pr0*dsqrt(1.d0+calph_p*ri(i,j,k))
     &       +vdts(i+1,j,k)*pr0*dsqrt(1.d0+calph_p*ri(i+1,j,k)))
     &    +anhf(j)*(vdts(i,j+1,k)*pr0*dsqrt(1.d0+calph_p*ri(i,j+1,k))
     &       +vdts(i+1,j+1,k)*pr0*dsqrt(1.d0+calph_p*ri(i+1,j+1,k))))
#endif
         vduv(i,j,k)=dmin1(vduv(i,j,k),vduvmx)
         vduv(i,j,k)=dmax1(vduv(i,j,k),vduvmn)*ex(i,j,k)
      enddo
      enddo
      enddo

      return
      end

