ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     sea-ice model
c
c     $Id: ice.F90 2 2007-11-15 13:07:08Z ishikawa $
c
c     ice_pre is called by surf-bulk.F
c     ice_fwd is called by tracli.F
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ice_pre

      use param
      use mod_mpi
      implicit none

#include "common.h"

      real*8 uicea(im,jm),uiceb(im,jm),aiceu(im,jm),
     &     vicea(im,jm),viceb(im,jm),uicep,vicep,
     &     dticor(jm),dticor2(jm)
      integer i,j,n

      call exch_t2d_n1p(aice,volice,1)
      call exch_t2d_s1p(aice,volice,2)
c
      call wait_t2d_n1p(aice,volice,1)
      call wait_t2d_s1p(aice,volice,2)
c
      call exch_t2d_e1p(aice,volice,3)
      call exch_t2d_w1p(aice,volice,4)
c
ccc   initialize
      do j = 1,jm
      do i = 1,im
         hiceu(i,j)=0.d0
      enddo
      enddo
c
      call wait_t2d_e1p(aice,volice,3)
c
      do j = 3,jml-2
      do i = 3,iml-2
!   hice of uv grid
         aiceu(i,j)=ex(i,j,1)*(ashf(j)*(aice(i,j)+aice(i+1,j))
     &        +anhf(j)*(aice(i,j+1)+aice(i+1,j+1)))/areauu(j)
c
         if(aice(i,j)*tex(i,j,1).gt.0.d0) then
            hiceu(i,j)=hiceu(i,j)+ashf(j)/areauu(j)
     &           *volice(i,j)/(areat(i,j,1)*aice(i,j))
         endif
         if(aice(i+1,j)*tex(i+1,j,1).gt.0.d0) then
            hiceu(i,j)=hiceu(i,j)+ashf(j)/areauu(j)
     &           *volice(i+1,j)/(areat(i+1,j,1)*aice(i+1,j))
         endif
         if(aice(i,j+1)*tex(i,j+1,1).gt.0.d0) then
            hiceu(i,j)=hiceu(i,j)+anhf(j)/areauu(j)
     &           *volice(i,j+1)/(areat(i,j+1,1)*aice(i,j+1))
         endif
         if(aice(i+1,j+1)*tex(i+1,j+1,1).gt.0.d0) then
            hiceu(i,j)=hiceu(i,j)+anhf(j)/areauu(j)
     &           *volice(i+1,j+1)/(areat(i+1,j+1,1)*aice(i+1,j+1))
         endif
c
         if(hiceu(i,j).lt.hmin)then
            hiceu(i,j) = 0.d0
            aiceu(i,j) = 0.d0
         endif
c
ccc   ice pressure for ice velocity prediction
c         pice(i,j)=tex(i,j,1)*5.d0*rho_i*1.e-3*grav*volice(i,j)
c     &        *exp(-20.d0*(1.d0-aice(i,j)))
c     &        /(areat(i,j,1)+1.d0-tex(i,j,1))
      enddo
      enddo
c
      do j=1,jml
         dticor(j)=.5d0*cor(j)*dtuice
         dticor2(j)=1.d0+dticor(j)*dticor(j)
      enddo
c
      do j = 1,jml
      do i = 1,iml
         if(hiceu(i,j).lt.hmin)then
            uice(i,j)=0.d0
            vice(i,j)=0.d0
         endif
c
         uiceb(i,j) = uice(i,j)
         viceb(i,j) = vice(i,j)
c
         uicea(i,j) = 0.d0
         vicea(i,j) = 0.d0
      enddo
      enddo
c
      do n=1,nint(dtts/dtuice)
c
      do j = 3,jml-2
      do i = 3,iml-2
      if(hiceu(i,j)*ex(i,j,1).gt.0.d0)then
         uicep=cd_a_i/cd_a_o*wsx(i,j)/(rho_i*1.d-3)/hiceu(i,j)
     &      +cor(j)*viceb(i,j)-dump_i*uiceb(i,j)
     &      -dsqrt((uiceb(i,j)-u(i,j,1))**2+(viceb(i,j)-v(i,j,1))**2)
     &      *(uiceb(i,j)-u(i,j,1))*rho_w*cd_i_o/rho_i/hiceu(i,j)
         vicep=cd_a_i/cd_a_o*wsy(i,j)/(rho_i*1.d-3)/hiceu(i,j)
     &      -cor(j)*uiceb(i,j)-dump_i*viceb(i,j)
     &      -dsqrt((uiceb(i,j)-u(i,j,1))**2+(viceb(i,j)-v(i,j,1))**2)
     &      *(viceb(i,j)-v(i,j,1))*rho_w*cd_i_o/rho_i/hiceu(i,j)
c
         uicea(i,j)=uiceb(i,j)
     &        +dtuice/dticor2(j)*(uicep+dticor(j)*vicep)
         vicea(i,j)=viceb(i,j)
     &        +dtuice/dticor2(j)*(vicep-dticor(j)*uicep)
      endif
      enddo
      enddo
c
      do j = 1,jml
      do i = 1,iml
         uiceb(i,j) = uicea(i,j)
         viceb(i,j) = vicea(i,j)
      enddo
      enddo
c
      enddo
c
      do j = 1,jml
      do i = 1,iml
      if(hiceu(i,j).ge.hmin)then
         uice(i,j) = uiceb(i,j)*ex(i,j,1)
         vice(i,j) = viceb(i,j)*ex(i,j,1)
      else
         uice(i,j)=0.d0
         vice(i,j)=0.d0
      endif
      enddo
      enddo
#ifdef NESTED
!    side boundary data for nest step
!  west
      if(ip_x .eq. 0) then
      do j=1,jml
      do i=1,2
         if(hiceu(i+2,j).ge.hmin) then
            uice(i+2,j)=uicewbt(i,j)*ex(i+2,j,1)
            vice(i+2,j)=vicewbt(i,j)*ex(i+2,j,1)
         else
            uice(i+2,j)=0.d0
            vice(i+2,j)=0.d0
         endif
      enddo
      enddo
      do j=1,jml
      do i=1,2
         uicewbt(i,j)=uicewbc(i,j,mb1)
     &           +(uicewbc(i,j,mb2)-uicewbc(i,j,mb1))*cbf
         vicewbt(i,j)=vicewbc(i,j,mb1)
     &           +(vicewbc(i,j,mb2)-vicewbc(i,j,mb1))*cbf
      enddo
      enddo
      endif
!east
      if(ip_x .eq. ipe-1) then
      do j=1,jml
      do i=1,2
         if(hiceu(iml-2-i,j).ge.hmin) then
            uice(iml-2-i,j)=uiceebt(i,j)*ex(iml-2-i,j,1)
            vice(iml-2-i,j)=viceebt(i,j)*ex(iml-2-i,j,1)
         else
            uice(iml-2-i,j)=0.d0
            vice(iml-2-i,j)=0.d0
         endif           
      enddo
      enddo
      do j=1,jml
      do i=1,2
         uiceebt(i,j)=uiceebc(i,j,mb1)
     &           +(uiceebc(i,j,mb2)-uiceebc(i,j,mb1))*cbf
         viceebt(i,j)=viceebc(i,j,mb1)
     &           +(viceebc(i,j,mb2)-viceebc(i,j,mb1))*cbf
      enddo
      enddo      
      endif
! south
      if(ip_y .eq. 0) then
      do i=1,iml
      do j=1,2
         if (hiceu(i,j+2).ge.hmin) then
            uice(i,j+2)=uicesbt(i,j)*ex(i,j+2,1)
            vice(i,j+2)=vicesbt(i,j)*ex(i,j+2,1)
         else
            uice(i,j+2)=0.d0
            vice(i,j+2)=0.d0
         endif
      enddo
      enddo
      do i=1,iml
      do j=1,2
         uicesbt(i,j)=uicesbc(i,j,mb1)
     &           +(uicesbc(i,j,mb2)-uicesbc(i,j,mb1))*cbf
         vicesbt(i,j)=vicesbc(i,j,mb1)
     &           +(vicesbc(i,j,mb2)-vicesbc(i,j,mb1))*cbf
      enddo
      enddo      
      endif
! north
      if(ip_y .eq. jpe-1) then
      do i=1,iml
      do j=1,2
         if(hiceu(i,jml-2-j).ge.hmin) then
            uice(i,jml-2-j)=uicenbt(i,j)*ex(i,jml-2-j,1)
            vice(i,jml-2-j)=vicenbt(i,j)*ex(i,jml-2-j,1)
         else
            uice(i,jml-2-j)=0.d0
            vice(i,jml-2-j)=0.d0
         endif
      enddo
      enddo
      do i=1,iml
      do j=1,2
         uicenbt(i,j)=uicenbc(i,j,mb1)
     &           +(uicenbc(i,j,mb2)-uicenbc(i,j,mb1))*cbf
         vicenbt(i,j)=vicenbc(i,j,mb1)
     &           +(vicenbc(i,j,mb2)-vicenbc(i,j,mb1))*cbf
      enddo
      enddo      
      endif
#endif

c
ccc   change wind stress for ocean under ice
      do j = 3,jml-2
      do i = 3,iml-2
      if(aiceu(i,j)*ex(i,j,1).gt.0.d0)then
         wsx(i,j)=wsx(i,j)*(1.d0-aiceu(i,j))
     &        +ex(i,j,1)*rho_w*1.d-3*cd_i_o*dsqrt(
     &        (uice(i,j)-u(i,j,1))**2+(vice(i,j)-v(i,j,1))**2
     &        )*(uice(i,j)-u(i,j,1))*aiceu(i,j)
         wsy(i,j)=wsy(i,j)*(1.d0-aiceu(i,j))
     &        +ex(i,j,1)*rho_w*1.d-3*cd_i_o*dsqrt(
     &        (uice(i,j)-u(i,j,1))**2+(vice(i,j)-v(i,j,1))**2
     &        )*(vice(i,j)-v(i,j,1))*aiceu(i,j)
      endif
      enddo
      enddo
c
c
ccc   sea-ice pressure
      do j = 3,jml-2
      do i = 3,iml-2
         dpdx(i,j,0)=ex(i,j,1)*rho_i*1.d-3*grav*dx2r*csr(j)*ddmnar
     &        *(volice(i+1,j+1)/(areat(i+1,j+1,1)+1.-tex(i+1,j+1,1))
     &        +volice(i+1,j)/(areat(i+1,j,1)+1.-tex(i+1,j,1))
     &        -volice(i,j+1)/(areat(i,j+1,1)+1.-tex(i,j+1,1))
     &        -volice(i,j)/(areat(i,j,1)+1.-tex(i,j,1)))
         dpdy(i,j,0)=ex(i,j,1)*rho_i*1.d-3*grav*dy2r*ddmnar
     &        *(volice(i+1,j+1)/(areat(i+1,j+1,1)+1.-tex(i+1,j+1,1))
     &        -volice(i+1,j)/(areat(i+1,j,1)+1.-tex(i+1,j,1))
     &        +volice(i,j+1)/(areat(i,j+1,1)+1.-tex(i,j+1,1))
     &        -volice(i,j)/(areat(i,j,1)+1.-tex(i,j,1)))
      enddo
      enddo
c
      call wait_t2d_w1p(aice,volice,4)
c
      return
      end
c
c
c------------------------------------------------------------------------
c
      subroutine ice_fwd

      use param
      use mod_mpi
      implicit none

#include "common.h"
      real*8,parameter :: aicemax=1.d0
      integer i,j,k
      real*8 wf_ice(im,jm),facw,face,uicew,uicee,vices,vicen,
     &     vol_di,t_frz,vol_oi,h_ice,vol_cr(im,jm),facup,facdn,
     &     wf1,wf2,tvt1,svt1,tvt2,svt2,hf_ice(im,jm),h_iceb(im,jm)
#ifdef NESTED
      integer jsn,jen
      real*8 adblc
#endif
      
      call exch_t2d_s1p(uice,vice,1)
      call wait_t2d_s1p(uice,vice,1)
      call exch_t2d_w1p(uice,vice,2)
c
      if(matsno.eq.1.and.matsn2.eq.0) then
         do j = 1,jml
         do i = 1,iml
            aiceb(i,j)=aice(i,j)
            voliceb(i,j)=volice(i,j)
         enddo
         enddo
      endif
c
      do j = 1,jml
      do i = 1,iml
         aicea(i,j)=0.d0
         volicea(i,j)=0.d0
         hf_ice(i,j) = 0.d0
         vol_cr(i,j) = 0.d0
         if(areat(i,j,1)*aiceb(i,j).ne.0.) then
            h_iceb(i,j) = voliceb(i,j)/(areat(i,j,1)*aiceb(i,j))
         else
            h_iceb(i,j) = 0.
         endif
         wf_ice(i,j)=0.
#ifdef DEBUG1
      if(nkai.ge.1670 .and. nkai .le. 1677 .and. ip.eq.14
     & .and. i.eq.16 .and. j.eq.25) then
        write(*,*) nkai,'init',aice(i,j),volice(i,j),h_iceb(i,j),
     &    aiceb(i,j),voliceb(i,j)
      endif
#endif
      enddo
      enddo
      


      call wait_t2d_w1p(uice,vice,2)
c
c
ccc   advection and diffusion of sea ice
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,1).eq.1.)then
c
      uicew=0.5d0*(uice(i-1,j-1)+uice(i-1,j))
      if(uicew.gt.0.d0) then
         facw =1.               !(1.d0+hupp)*0.5d0
         face =0.               !(1.d0-hupp)*0.5d0
         if(aiceb(i,j).gt.0.99d0)facw=0.
      else
         facw =0.               !(1.d0-hupp)*0.5d0
         face =1.               !(1.d0+hupp)*0.5d0
         if(aiceb(i-1,j).gt.0.99d0)face=0.
      endif
c
      volicea(i,j)=volicea(i,j)
     &     -hdice*(dxr*cstr(j))**2
     &     *(voliceb(i,j)/areat(i,j,1)
     &     -voliceb(i-1,j)/(areat(i-1,j,1)+1.-tex(i-1,j,1)))
     &     *(ex(i-1,j-1,1)*anhf(j-1)+ex(i-1,j,1)*ashf(j))*2.
     &     *tex(i-1,j,1)
      if(facw*voliceb(i-1,j)+face*voliceb(i,j).gt.0.)then
         volicea(i,j)=volicea(i,j)
     &        +uicew*(facw*voliceb(i-1,j)
     &        +face*voliceb(i,j))*dxr*cstr(j)
     &        *tex(i-1,j,1)
      endif
      uicee=0.5d0*(uice(i,j-1)+uice(i,j))
      if(uicee.gt.0.d0)then
         facw =1.               !(1.d0+hupp)*0.5d0
         face =0.               !(1.d0-hupp)*0.5d0
         if(aiceb(i+1,j).gt.0.99d0)facw=0.
      else
         facw =0.               !(1.d0-hupp)*0.5d0
         face =1.               !(1.d0+hupp)*0.5d0
         if(aiceb(i,j).gt.0.99d0)face=0.
      endif
      volicea(i,j)=volicea(i,j)
     &     +hdice*(dxr*cstr(j))**2
     &     *(voliceb(i+1,j)/(areat(i+1,j,1)+1.-tex(i+1,j,1))
     &     -voliceb(i,j)/areat(i,j,1))
     &     *(ex(i,j-1,1)*anhf(j-1)+ex(i,j,1)*ashf(j))*2.
     &     *tex(i+1,j,1)
      if(facw*voliceb(i,j)+face*voliceb(i+1,j).gt.0.)then
         volicea(i,j)=volicea(i,j)
     &        -uicee*(facw*voliceb(i,j)
     &        +face*voliceb(i+1,j))*dxr*cstr(j)
     &        *tex(i+1,j,1)
      endif

      vices=0.5d0*(vice(i-1,j-1)+vice(i,j-1))
      if(vices.gt.0.d0) then
         facw =1.               !(1.d0+hupp)*0.5d0
         face =0.               !(1.d0-hupp)*0.5d0
         if(aiceb(i,j).gt.0.99d0)facw=0.
      else
         facw =0.               !(1.d0-hupp)*0.5d0
         face =1.               !(1.d0+hupp)*0.5d0
         if(aiceb(i,j-1).gt.0.99d0)face=0.
      endif
c
      volicea(i,j)=volicea(i,j)
     &     -hdice*dyr*dyr
     &     *(voliceb(i,j)/areat(i,j,1)
     &     -voliceb(i,j-1)/(areat(i,j-1,1)+1.-tex(i,j-1,1)))
     &     *(ex(i-1,j-1,1)+ex(i,j-1,1))*(ashf(j-1)+anhf(j-1))
     &     *tex(i,j-1,1)
      if(facw*voliceb(i,j-1)+face*voliceb(i,j).gt.0.)then
         volicea(i,j)=volicea(i,j)
     &        +vices*(facw*voliceb(i,j-1)
     &        +face*voliceb(i,j))*dyr
     &        *tex(i,j-1,1)
      endif

      vicen=0.5d0*(vice(i-1,j)+vice(i,j))
      if(vicen.gt.0.d0)then
         facw =1.               !(1.d0+hupp)*0.5d0
         face =0.               !(1.d0-hupp)*0.5d0
         if(aiceb(i,j+1).gt.0.99d0)facw=0.
      else
         facw =0.               !(1.d0-hupp)*0.5d0
         face =1.               !(1.d0+hupp)*0.5d0
         if(aiceb(i,j).gt.0.99d0)face=0.
      endif
c
      volicea(i,j)=volicea(i,j)
     &     +hdice*dyr*dyr
     &     *(voliceb(i,j+1)/(areat(i,j+1,1)+1.-tex(i,j+1,1))
     &     -voliceb(i,j)/areat(i,j,1))
     &     *(ex(i-1,j,1)+ex(i,j,1))*(ashf(j)+anhf(j))
     &     *tex(i,j+1,1)
      if(facw*voliceb(i,j)+face*voliceb(i,j+1).gt.0.)then
         volicea(i,j)=volicea(i,j)
     &        -vicen*(facw*voliceb(i,j)+face*voliceb(i,j+1))*dyr
     &        *tex(i,j+1,1)
      endif


#ifdef NESTED
      adblc=1.d0
!   west
      if(ip_x.eq.0) then
         if(i.ge.5 .and. i.le.indg+4) then
            volicea(i,j)=volicea(i,j)
     $           +chfice*(volicewbt(i-2,j)-volice(i,j))
     $           /(areat(i,j,1)+1.d0-tex(i,j,1))
            adblc=0.d0
         endif
      endif

!   east
      if(ip_x.eq.ipe-1) then
         if(i.le.iml-4 .and. i.ge.iml-indg-3) then
            volicea(i,j)=volicea(i,j)
     $           +chfice*(voliceebt(iml-i-1,j)-volice(i,j))
     $           /(areat(i,j,1)+1.d0-tex(i,j,1))
            adblc=0.d0
         endif
      endif

!   south
      if(ip_y.eq.0) then
         if(j.ge.5 .and. j.le.jndg+4) then
            volicea(i,j)=volicea(i,j)
     $           +chfice*(volicesbt(i,j-2)-volice(i,j))*adblc
     $           /(areat(i,j,1)+1.d0-tex(i,j,1))
         endif
      endif

!   north
      if(ip_y.eq.jpe-1) then
         if(j.le.jml-4 .and. j.ge.jml-jndg-3) then
            volicea(i,j)=volicea(i,j)
     $           +chfice*(volicenbt(i,jml-j-1)-volice(i,j))*adblc
     $           /(areat(i,j,1)+1.d0-tex(i,j,1))
         endif
      endif
#endif

#ifdef DEBUG1
      if(nkai.ge.1670 .and. nkai .le. 1677 .and. ip.eq.14
     & .and. i.eq.16 .and. j.eq.25) then
        write(*,*) nkai,'nest',volicea(i,j)
      endif
#endif
ccc   dynamics of ice
c
      volicea(i,j)=tex(i,j,1)*(voliceb(i,j)+volicea(i,j)*c2dtts)
      volicea(i,j)=dmax1(volicea(i,j),0.d0)
      
      if(aiceb(i,j).eq.0. .and. volicea(i,j).ge.0.)then
         aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmin)
#ifdef DEBUG1
      if(nkai.ge.1670 .and. nkai .le. 1677 .and. ip.eq.14
     & .and. i.eq.16 .and. j.eq.25) then
        write(*,*) nkai,'dyn1',volicea(i,j),aicea(i,j),hmin,areat(i,j,1)
      endif
#endif
      else
         aicea(i,j)=aiceb(i,j)*volicea(i,j)/voliceb(i,j)
#ifdef DEBUG1
      if(nkai.ge.1670 .and. nkai .le. 1677 .and. ip.eq.14
     & .and. i.eq.16 .and. j.eq.25) then
        write(*,*) nkai,'dyn2',volicea(i,j),aicea(i,j),voliceb(i,j),
     &   aiceb(i,j)
      endif
#endif
      endif
c
ccc   growth and melting of sea ice
c
      t_frz=sa(i,j,1)*(a1_frz+a2_frz*dsqrt(sa(i,j,1))
     &   +a3_frz*sa(i,j,1))
c
ccc   direct flux
      vol_di=-areat(i,j,1)*netq_i(i,j)*c2dtts/(rho_i*c_p_i*1.d-6)
      vol_di=dmax1(vol_di,-volicea(i,j))
      hf_ice(i,j) = -t_frz*vol_di*rho_i/rho_w
c
ccc   through the ocean
c
      if(aicea(i,j).gt.aicemin.or.ta(i,j,1).lt.t_frz) then
         vol_oi=volt(i,j,1)*c_p_w*rho_w*(t_frz-ta(i,j,1))
     &        /(c_p_i*rho_i)
         vol_oi=dmax1(vol_oi,-volicea(i,j))
c         vol_oi=dmin1(vol_oi,areat(i,j,1)*hmax-volicea(i,j))
         hf_ice(i,j) = hf_ice(i,j)
     &        +vol_oi*rho_i*c_p_i/rho_w/c_p_w
      else
         vol_oi = 0.
      endif
c
      if(volicea(i,j).eq.0.d0)then
         aicea(i,j)=vol_oi/(areat(i,j,1)*hmin+1.-tex(i,j,1))
         t_top(i,j) = t_frz
      else
         aicea(i,j)=aicea(i,j)*(volicea(i,j)+vol_oi)/volicea(i,j)
      endif
      volicea(i,j)=tex(i,j,1)*(volicea(i,j)+vol_oi+vol_di)
c
ccc   constraint for ice
c
      if(aicea(i,j).lt.aicemin.or.volicea(i,j).le.0.d0) then
ccc   all ice melt
         vol_cr(i,j) = -volicea(i,j)
         aicea(i,j) = 0.
         volicea(i,j) = 0.
c
      else
ccc   area over
         if(aicea(i,j) > aicemax) then
            aicea(i,j) = aicemax
         endif
c
ccc   ice thickness condition
         h_ice=volicea(i,j)/(aicea(i,j)*areat(i,j,1)+1.-tex(i,j,1))
     &        *tex(i,j,1)
ccc   thin ice
         if(h_ice.lt.hmin) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmin+1.-tex(i,j,1))
            if(aicea(i,j).lt.aicemin) then
               vol_cr(i,j) = -volicea(i,j)
               aicea(i,j) = 0.
               volicea(i,j) = 0.
            endif
         endif
c
ccc   thick ice
         if(h_ice.gt.hmax) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmax+1.-tex(i,j,1))
            if(aicea(i,j).gt.1.d0) then
               vol_cr(i,j) = areat(i,j,1)*hmax-volicea(i,j)
               aicea(i,j) = 1.d0
               volicea(i,j) = areat(i,j,1)*hmax
            endif
         endif
c
      endif
c
      hf_ice(i,j)=hf_ice(i,j)+vol_cr(i,j)*rho_i*c_p_i/rho_w/c_p_w
c
ccc   fresh water flux
      wf_ice(i,j)=-(vol_di+vol_oi+vol_cr(i,j))*rho_i/rho_w
      hcla(i,j) = hcla(i,j) + wf_ice(i,j)/(areat(i,j,1)+1.-tex(i,j,1))
c
ccc   sigma coordinate
c      if(wf_ice(i,j).gt.0.d0) then
      if(wf_ice(i,j).lt.0.d0) then
         facup = (1.d0-vupp)*0.5d0
         facdn = (1.d0+vupp)*0.5d0
      else
         facup = (1.d0+vupp)*0.5d0
         facdn = (1.d0-vupp)*0.5d0
      endif
c
      wf1 = -wf_ice(i,j)*(1.d0-dep(2)/depadp)
      tvt1=wf1*(facup*ta(i,j,1)+facdn*ta(i,j,2))
      svt1=wf1*(facup*sa(i,j,1)+facdn*sa(i,j,2))
c
ccc   salt flux
      sa(i,j,1) = tex(i,j,1)*(sa(i,j,1)*volt(i,j,1)
     &     +s_ice*wf_ice(i,j)*rho_w*1.d-3+svt1)
     &   /(volt(i,j,1)+wf_ice(i,j)*dzt(i,j,1)/depadp+1.-tex(i,j,1))
c
ccc   heat flux
      ta(i,j,1) = tex(i,j,1)*(ta(i,j,1)*volt(i,j,1)+hf_ice(i,j)+tvt1)
     &   /(volt(i,j,1)+wf_ice(i,j)*dzt(i,j,1)/depadp+1.-tex(i,j,1))
      endif
      enddo
      enddo
      
#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do j=3,jml-2
      do i=1,indg+2
         volicewbt(i,j)=(volicewbc(i,j,mb1)
     $        +(volicewbc(i,j,mb2)-volicewbc(i,j,mb1))*cbf)
      enddo
      do i=3,4
         volicea(i,j)=volicewbt(i-2,j)*tex(i,j,1)
         aicea(i,j)=(aicewbc(i-2,j,mb1)
     $        +(aicewbc(i-2,j,mb2)-aicewbc(i-2,j,mb1))*cbf)*tex(i,j,1)
c
         if(aicea(i,j).lt.aicemin.or.volicea(i,j).le.0.d0) then
            aicea(i,j) = 0.
            volicea(i,j) = 0.
         else
            if(aicea(i,j)>aicemax) then
               aicea(i,j) = aicemax
            endif
            h_ice=volicea(i,j)/(aicea(i,j)*areat(i,j,1)+1.-tex(i,j,1))
     &           *tex(i,j,1)
            if(h_ice.lt.hmin) then
               aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmin+1.-tex(i,j,1))
               if(aicea(i,j).lt.aicemin) then
                  aicea(i,j) = 0.
                  volicea(i,j) = 0.
               endif
            endif
            if(h_ice.gt.hmax) then
               aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmax+1.-tex(i,j,1))
               if(aicea(i,j).gt.1.d0) then
                  aicea(i,j) = 1.d0
                  volicea(i,j) = areat(i,j,1)*hmax
               endif
            endif
         endif
c
      enddo
      enddo      
      endif

! east
      if( ip_x .eq. ipe-1) then
      do j=3,jml-2
      do i=1,indg
         voliceebt(i,j)=voliceebc(i,j,mb1)
     &        +(voliceebc(i,j,mb2)-voliceebc(i,j,mb1))*cbf
      enddo
      do i=iml-3,iml-2
         volicea(i,j)=voliceebt(iml-i-1,j)*tex(i,j,1)
         aicea(i,j)=(aiceebc(iml-1-i,j,mb1)
     &    +(aiceebc(iml-1-i,j,mb2)-aiceebc(iml-1-i,j,mb1))*cbf)
     $           *tex(i,j,1)

         if(aicea(i,j).lt.aicemin.or.volicea(i,j).le.0.d0) then
            aicea(i,j) = 0.
            volicea(i,j) = 0.
         else
            if(aicea(i,j) > aicemax) then
               aicea(i,j) = aicemax
            endif
            h_ice=volicea(i,j)/(aicea(i,j)*areat(i,j,1)+1.-tex(i,j,1))
     &           *tex(i,j,1)
         if(h_ice.lt.hmin) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmin+1.-tex(i,j,1))
            if(aicea(i,j).lt.aicemin) then
               aicea(i,j) = 0.
               volicea(i,j) = 0.
            endif
         endif
         if(h_ice.gt.hmax) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmax+1.-tex(i,j,1))
            if(aicea(i,j)>aicemax) then
               aicea(i,j) = aicemax
               volicea(i,j) = areat(i,j,1)*hmax*aicemax
            endif
         endif
      endif
c
      enddo
      enddo      
      endif

! south
      if(ip_y .eq. 0) then
      do i=1,iml
      do j=1,indg+2
          volicesbt(i,j)=volicesbc(i,j,mb1)
     $        +(volicesbc(i,j,mb2)-volicesbc(i,j,mb1))*cbf
      enddo
      do j=3,4
         volicea(i,j)=volicesbt(i,j-2)*tex(i,j,1)
         aicea(i,j)=(aicesbc(i,j-2,mb1)
     $        +(aicesbc(i,j-2,mb2)-aicesbc(i,j-2,mb1))*cbf)*tex(i,j,1)
c
         if(aicea(i,j).lt.aicemin.or.volicea(i,j).le.0.d0) then
            aicea(i,j) = 0.
            volicea(i,j) = 0.
         else
            if(aicea(i,j).gt.1.d0) then
               aicea(i,j) = 1.d0
            endif
            h_ice=volicea(i,j)/(aicea(i,j)*areat(i,j,1)+1.-tex(i,j,1))
     &           *tex(i,j,1)
         if(h_ice.lt.hmin) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmin+1.-tex(i,j,1))
            if(aicea(i,j).lt.aicemin) then
               aicea(i,j) = 0.
               volicea(i,j) = 0.
            endif
         endif
         if(h_ice.gt.hmax) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmax+1.-tex(i,j,1))
            if(aicea(i,j).gt.1.d0) then
               aicea(i,j) = 1.d0
               volicea(i,j) = areat(i,j,1)*hmax
            endif
         endif
      endif
c
      enddo
      enddo      
      endif

! north
      if(ip_y .eq. jpe-1) then      
      do i=1,iml
      do j=1,jndg+2
         volicenbt(i,j)=volicenbc(i,j,mb1)
     $        +(volicenbc(i,j,mb2)-volicenbc(i,j,mb1))*cbf
      enddo
      do j=jml-3,jml-2
         volicea(i,j)=volicenbt(i,jml-1-j)*tex(i,j,1)
         aicea(i,j)=(aicenbc(i,jml-1-j,mb1)
     $    +(aicenbc(i,jml-1-j,mb2)-aicenbc(i,jml-1-j,mb1))*cbf)
     $        *tex(i,j,1)
c
         if(aicea(i,j).lt.aicemin.or.volicea(i,j).le.0.d0) then
            aicea(i,j) = 0.
            volicea(i,j) = 0.
         else
            if(aicea(i,j).gt.1.d0) then
               aicea(i,j) = 1.d0
            endif
            h_ice=volicea(i,j)/(aicea(i,j)*areat(i,j,1)+1.-tex(i,j,1))
     &           *tex(i,j,1)
         if(h_ice.lt.hmin) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmin+1.-tex(i,j,1))
            if(aicea(i,j).lt.aicemin) then
               aicea(i,j) = 0.
               volicea(i,j) = 0.
            endif
         endif
         if(h_ice.gt.hmax) then
            aicea(i,j)=volicea(i,j)/(areat(i,j,1)*hmax+1.-tex(i,j,1))
            if(aicea(i,j).gt.1.d0) then
               aicea(i,j) = 1.d0
               volicea(i,j) = areat(i,j,1)*hmax
            endif
         endif
      endif

      enddo
      enddo
      endif
#endif
c
ccc   vertical advection due to the sigma coordiate
c
      do j = 3,jml-2
      do k = 2,kadp
      do i = 3,iml-2
         if(wf_ice(i,j).gt.0.d0) then
            facup = (1.d0-vupp)*0.5d0
            facdn = (1.d0+vupp)*0.5d0
         else
            facup = (1.d0+vupp)*0.5d0
            facdn = (1.d0-vupp)*0.5d0
         endif
c
         wf1 = -wf_ice(i,j)*(1.d0-dep(k)/depadp)
         tvt1=wf1*(facup*ta(i,j,k-1)+facdn*ta(i,j,k))
         svt1=wf1*(facup*sa(i,j,k-1)+facdn*sa(i,j,k))
c
         wf2 = -wf_ice(i,j)*(1.d0-dep(k+1)/depadp)
         tvt2=wf2*(facup*ta(i,j,k)+facdn*ta(i,j,k+1))
         svt2=wf2*(facup*sa(i,j,k)+facdn*sa(i,j,k+1))
c
         ta(i,j,k) = (ta(i,j,k)*volt(i,j,k)+tvt2-tvt1)
     &     /(volt(i,j,k)+wf_ice(i,j)*dzt(i,j,k)/depadp+1.-tex(i,j,k))
         sa(i,j,k) = (sa(i,j,k)*volt(i,j,k)+svt2-svt1)
     &     /(volt(i,j,k)+wf_ice(i,j)*dzt(i,j,k)/depadp+1.-tex(i,j,k))
c
#ifdef DEBUG1
         if(sa(i,j,k)+1.-tex(i,j,k).le.0.d0) then
            write(*,*) 'S# (ice) ',nkai,lw(ip_x)-1+i,ls(ip_y)+j-1,k,
     &           sa(i,j,k),tex(i,j,k),' -> stop; local ij',i,j
!            stop
         endif
#endif
      enddo
      enddo
      enddo
c
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,k).eq.1.d0)then
         t_frz=s(i,j,k)*(a1_frz+a2_frz*dsqrt(s(i,j,k))
     &        +a3_frz*s(i,j,k))*tex(i,j,k)
         ta(i,j,k)=dmax1(ta(i,j,k),t_frz)*tex(i,j,k)
      endif
      enddo
      enddo
      enddo
c
#ifdef DEBUG1
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,k).ne.0.)then
         if(ta(i,j,k).ge.35..or.ta(i,j,k).le.-3.0) then
            write(*,*)'T# (ice),nkai= ',nkai,'ip,ip_x,ip_y',
     &      ip,ip_x,ip_y,'i,j',i,j,'ig,jg',
     &           lw(ip_x)+i-1,ls(ip_y)+j-1,'k ',k,
     &           'ta/t/tb ',ta(i,j,k),t(i,j,k),tb(i,j,k)
            write(*,*)'  i-1/+1/j-1/+1 ',
     &           t(i-1,j,k),t(i+1,j,k),t(i,j-1,k),t(i,j+1,k)
            write(*,*)'  k-1/+1/r/g    ',
     &           t(i,j,k-1),t(i,j,k+1),tref(i,j,k),gref(i,j,k)
            write(*,*)'  w/vddt ',wl(i,j,k)/areat(i,j,k),
     &           wl(i,j,k+1)/(areat(i,j,k+1)+1.-tex(i,j,k+1)),
     &           vddt(i,j,k),vddt(i,j,k+1)
         endif
         if(sa(i,j,k).ge.40..or.sa(i,j,k).le.5.) then
            write(*,*) 'S# (ice)nkai= ',nkai,'ip,ip_x,ip_y',
     &      ip,ip_x,ip_y,'i,j',i,j,'ig,jg',
     &           lw(ip_x)+i-1,ls(ip_y)+j-1,'k ',k,
     &           'sa/s/sb ',sa(i,j,k),s(i,j,k),sb(i,j,k)
            write(*,*)'  i-1/+1/j-1/+1 ',
     &           s(i-1,j,k),s(i+1,j,k),s(i,j-1,k),s(i,j+1,k)
            write(*,*)'  k-1/+1/r/g    ',
     &           s(i,j,k-1),s(i,j,k+1),sref(i,j,k),gref(i,j,k)
            write(*,*)'  w/vdds ',wl(i,j,k)/areat(i,j,k),
     &           wl(i,j,k+1)/(areat(i,j,k+1)+1.-tex(i,j,k+1)),
     &           vdds(i,j,k),vdds(i,j,k+1)
         endif
      endif
      enddo
      enddo
      enddo
#endif
#ifdef DEBUG1
      if(nkai.ge.1670 .and. nkai .le. 1677 .and. ip.eq.14) then
        i = 16
        j = 25
        write(*,*) nkai,'last',volicea(i,j),aicea(i,j)
      endif
#endif
c
ccc   transport things to prepare next step
c
      if(matsno.eq.1.or.matsn2.eq.1) then
         do j = 1,jml
         do i = 1,iml
            aice(i,j)=tex(i,j,1)*aicea(i,j)
            volice(i,j)=tex(i,j,1)*volicea(i,j)
         enddo
         enddo
      else
         do j = 1,jml
         do i = 1,iml 
            aiceb(i,j)=aice(i,j)
            voliceb(i,j)=volice(i,j)
            aice(i,j)=tex(i,j,1)*aicea(i,j)
            volice(i,j)=tex(i,j,1)*volicea(i,j)
         enddo
         enddo
      endif

      return
      end

