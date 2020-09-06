ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine tracli                                             c
c
c  $Id: tracli.F90 13 2008-12-12 10:53:49Z ishikawa $
c                                                                   c
c     This is the first part of the main computation in each time   c
c     step. time changes of tracers (temperature and salinity) due  c
c     to advection and diffusion, and those in the baclinic part of c
c     velocity (deviations from vertical-mean) are calculated. the  c
c     generalized arakawa scheme is used for momentum advection.    c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tracli

      use param
      use mod_mpi
#ifdef BNDTIDE
      use mod_time
      use mod_nao
#endif
      implicit none

#include "common.h"
c
      real*8 sroot,botvel,dpth,swpn,texex,
     &     fcgini,factws,fclfr,t_frz
      real*8 tee,see,tww,sww,tnn,snn,tss,sss,cue,cve,cun,cvn,
     &  tc20,tc02,sc20,sc02,tc11,sc11,tc01,sc01,
     &  td20,td02,sd02,sd20,td11,sd11,td10,sd10,
     &  tc00,tc10,td00,td01,sc00,sc10,sd00,sd01
      real*8 tuu,tdd,suu,sdd,tc2,sc2,tc0,sc0,tc1,sc1,wc0,
     & tvt1,svt1,facup,facdn
      integer i,j,k

      real*8 fact1,wuv,uvis_l,vvis_l,uvis_b,vvis_b,uadv,vadv,uvt,vvt
      real*8 defrate_sw,defrate_nw,defrate_se,defrate_ne
      real*8 stbm,stbdf,stbdh,stbuf,stbuh,stb0,stb1,stb2
      real*8 gref1
c
#ifdef NESTED
      real*8 factu,facth
      integer jsn,jen
#endif

#ifdef BNDTIDE
      real*8 h_tide(im,jm),gae
      type(Tdate) :: date_tide
#endif

!   for implicit vertical diffusion
      real*8 tria1(im,jm,km),trib1(im,jm,km),tric1(im,jm,km),
     &  tribet1(im,jm),tria2(im,jm,km),
     &  trib2(im,jm,km),tric2(im,jm,km),tribet2(im,jm),
     &  trigam2(im,jm,km),
     &  tria1t(im,jm,km),trib1t(im,jm,km),tric1t(im,jm,km),
     &  tria1s(im,jm,km),trib1s(im,jm,km),tric1s(im,jm,km),
     &  tribet1t(im,jm),trigam1t(im,jm,km),
     &  tribet1s(im,jm),trigam1s(im,jm,km)

      call exch_t2d_n1(wflux,25)
      call exch_t2d_s1(wflux,36)

      call exch_t3d_n1p(t,s,2,3)
      call exch_t3d_s1p(t,s,12,13)

      call exch_t3d_n1p(u,v,8,9)
      call exch_t3d_s1p(u,v,10,11)

      call exch_t2d_n1p(sfund,sfvnd,6)
      call exch_t2d_s1p(sfund,sfvnd,7)

#ifdef BISMFRIC
      call exch_t2d_n2p(sfund,sfvnd,29)
      call exch_t2d_s2p(sfund,sfvnd,30)
#endif

      call wait_t2d_n1(wflux,25)
      call wait_t2d_s1(wflux,36)
      call exch_t2d_w1(wflux,43)
      call exch_t2d_e1(wflux,44)

      call wait_t3d_n1p(t,s,2,3)
      call wait_t3d_s1p(t,s,12,13)
      call exch_t3d_e1p(t,s,45,46)
      call exch_t3d_w1p(t,s,47,48)

      call wait_t3d_n1p(u,v,8,9)
      call wait_t3d_s1p(u,v,10,11)
      call exch_t3d_e1p(u,v,78,79)
      call exch_t3d_w1p(u,v,56,57)

      call wait_t2d_n1p(sfund,sfvnd,6)
      call wait_t2d_s1p(sfund,sfvnd,7)

#ifdef BISMFRIC
      call wait_t2d_n2p(sfund,sfvnd,29)
      call wait_t2d_s2p(sfund,sfvnd,30)
      call exch_t2d_e2p(sfund,sfvnd,64)
      call exch_t2d_w2p(sfund,sfvnd,65)
#endif

      call exch_t2d_e1p(sfund,sfvnd,58)
      call exch_t2d_w1p(sfund,sfvnd,59)


!   matsuno scheme 1st step

      if(matsno.eq.1.and.matsn2.eq.0) then
         do k = 1,km
         do j = 1,jml
         do i = 1,iml
            ub(i,j,k)=u(i,j,k)
            vb(i,j,k)=v(i,j,k)
            tb(i,j,k)=t(i,j,k)
            sb(i,j,k)=s(i,j,k)
         enddo
         enddo
         enddo

         do j = 1,jml
         do i = 1,im
            hclb(i,j) = hcl(i,j)
         enddo
         enddo
      endif

      call exch_t2d_n1p(hcl,hclb,14)
      call exch_t2d_s1p(hcl,hclb,15)

      call exch_t3d_n2p(tb,sb,26,27)
      call exch_t3d_n1p(tb,sb,16,17)
      call exch_t3d_s1p(tb,sb,18,19)

      call exch_t3d_n1p(ub,vb,20,21)
      call exch_t3d_s1p(ub,vb,22,23)

#ifdef BISMFRIC
      call exch_t3d_n2p(ub,vb,31,32)
      call exch_t3d_s2p(ub,vb,33,34)
#endif

      call wait_t2d_n1p(hcl,hclb,14)
      call wait_t2d_s1p(hcl,hclb,15)
      call exch_t2d_w1p(hcl,hclb,41)
      call exch_t2d_e1p(hcl,hclb,42)

      call wait_t3d_n2p(tb,sb,26,27)
      call wait_t3d_n1p(tb,sb,16,17)
      call wait_t3d_s1p(tb,sb,18,19)
      call exch_t3d_e2p(tb,sb,49,50)
      call exch_t3d_e1p(tb,sb,51,52)
      call exch_t3d_w1p(tb,sb,53,54)

      call wait_t3d_n1p(ub,vb,20,21)
      call wait_t3d_s1p(ub,vb,22,23)

#ifdef BISMFRIC
      call wait_t3d_n2p(ub,vb,31,32)
      call wait_t3d_s2p(ub,vb,33,34)
      call exch_t3d_e2p(ub,vb,66,67)
      call exch_t3d_w2p(ub,vb,68,69)
#endif

      call exch_t3d_e1p(ub,vb,60,61)
      call exch_t3d_w1p(ub,vb,62,63)





!   initialize

      do k = 0,km+1
      do j = 1,jm
      do i = 1,im
         ua(i,j,k)=0.d0
         va(i,j,k)=0.d0
         ta(i,j,k)=0.d0
         sa(i,j,k)=0.d0

         ustar(i,j,k) = 0.d0
         vstar(i,j,k) = 0.d0
         dpdx(i,j,k) = 0.d0
         dpdy(i,j,k) = 0.d0

         tust(i,j,k) = 0.d0
         tvst(i,j,k) = 0.d0
         sust(i,j,k) = 0.d0
         svst(i,j,k) = 0.d0
      enddo
      enddo
      enddo

      do j = 1,jm
      do i = 1,im
         zu(i,j)=0.d0
         zv(i,j)=0.d0
         sfu(i,j) = 0.d0
         sfv(i,j) = 0.d0

         ta1toice(i,j)=0.d0
      enddo
      enddo

#ifdef NESTED
#ifdef REST_INI
#ifdef ROKKA
      factu = dmin1(dble(nkai)/dble(nxmonth),1.d0)
#endif
#ifdef NWNPAC
c      factu = dmin1(dble(nkai)/dble(3*nxyear),1.d0)
      factu = dmin1(dble(nkai)/dble(nxmonth),1.d0)
c      factu = dmin1(dble(nkai)/dble(3*nxday),1.d0)
#endif
#ifdef JP68M
      factu = dmin1(dble(nkai)/dble(nxmonth),1.d0)
#endif
#else
      factu=1.d0
#endif
#endif

#ifdef BNDTIDE
      date_tide = date_now + diff_dtuv
      call waveheight(date_tide,h_tide)
!     write(*,*) 'tracli :',ip,maxval(h_tide),minval(h_tide)
#endif !BNDTIDE



!   facters for calm spin-up
!   factws: wind stress
!   fcgini: initial restoring
!   fclfr : laplacian friction
      factws=1.d0
      fcgini=0.d0
#ifdef BISMFRIC
      fclfr=0.d0
#else
      fclfr=1.d0
#endif
#ifdef PC68M
       factws = dmin1(
     &  dble(nsec-181*24*3600)/dble( (365-181)*24*3600 ),1.d0)
!       factws = 0.
!      factws=dmin1(dble(nkai)/dble(nxyear),1.d0)
!      factws=dmin1(dble(nkai)/dble(nxmonth),1.d0)
!      fcgini=dmax1( dble(nxyear-nkai)/dble(nxyear),0.d0 )
	  fcgini=0.
c      fcgini=dmax1( dble(nxmonth-nkai)/dble(nxmonth),0.d0 )
!      fclfr=dmax1( dble(nxmonth-nkai)/dble(nxmonth),0.d0 )
#endif
#ifdef ROKKA
      factws=factu
!      fcgini=dmax1( dble(nxyear-nkai)/dble(nxyear),0.d0 )
!      fclfr=dmax1( dble(nxmonth-nkai)/dble(nxmonth),0.d0 )
      fcgini=0.d0
#ifdef MASKLVIS
      fclfr=1.d0
#else
      fclfr=0.d0
#endif
#endif
#ifdef BNDTIDE
       fclfr = 1.d0
#endif
#ifdef NWNPAC
      factws=factu
      fcgini=dmax1( dble(nxday*10-nkai)/dble(nxday*10),0.d0 )
      fclfr=0.d0
!      fclfr=dmax1( dble(nxmonth-nkai)/dble(nxmonth),0.d0 )
#endif
#ifdef JP68M
      factws=factu
      fcgini=1.d0
      fclfr=1.d0
!      factws=factu
!      fcgini=dmax1( dble(nxyear-nkai)/dble(nxyear),0.d0 )
!      fclfr=dmax1( dble(nxmonth-nkai)/dble(nxmonth),0.d0 )
#endif


!   next step hcl including water flux

      call wait_t2d_w1p(hcl,hclb,41)
      call wait_t2d_e1p(hcl,hclb,42)

      call wait_t2d_w1(wflux,43)
      call wait_t2d_e1(wflux,44)


#ifdef BNDTIDE
       wflux(:,:) = 0.d0
#endif
      do j = 2,jml-1
      do i = 2,iml-1
         hcla(i,j)=hclb(i,j)+tex(i,j,1)*c2dtts
     &   *( wflux(i,j)+wl(i,j,1)/(areat(i,j,1)+1.d0-tex(i,j,1)) )
      enddo
      enddo
#ifdef NESTED


! south
      if(ip_y .eq. 0) then
      do i=1,iml
      do j=1,jndg+2
         hclasbt(i,j)=hclsbc(i,j,mb1)
     $        +(hclsbc(i,j,mb2)-hclsbc(i,j,mb1))*cbf
      enddo
      enddo

      do i=3,iml-2
      do j=3,jndg+2
#ifdef BNDTIDE
       hclsbt(i,j) = h_tide(i,j+2)
#endif
         hcla(i,j+2)=(hcla(i,j+2)
     &         + chfhcl*c2dtts*(hclsbt(i,j)-hcl(i,j+2)))*tex(i,j+2,1)
         ! hcla(i,j+2)=hclsbt(i,j)*tex(i,j+2,1)
      enddo
      enddo

      do i=1,iml
      do j=1,2
#ifdef BNDTIDE
       hclasbt(i,j) = h_tide(i,j+2)
#endif
         hcla(i,j+2)=hclasbt(i,j)*tex(i,j+2,1)
      enddo
      enddo
! avoid double counting
        jsn = jndg+5
      else
        jsn = 3
      endif


!north
      if(ip_y .eq. jpe-1) then
      do i=1,iml
      do j=1,jndg+2
         hclanbt(i,j)=hclnbc(i,j,mb1)
     &        +(hclnbc(i,j,mb2)-hclnbc(i,j,mb1))*cbf
      enddo
      enddo

      do i=3,iml-2
      do j=3,jndg+2
#ifdef BNDTIDE
         hclnbt(i,j) =  h_tide(i,jml-1-j)
#endif
         hcla(i,jml-1-j)=(hcla(i,jml-1-j)
     &  + chfhcl*c2dtts*(hclnbt(i,j)-hcl(i,jml-1-j)))*tex(i,jml-1-j,1)
         ! hcla(i,jml-1-j)=hclnbt(i,j)*tex(i,jml-1-j,1)
      enddo
      enddo

      do i=1,iml
      do j=1,2
#ifdef BNDTIDE
         hclanbt(i,j) = h_tide(i,jml-1-j)
#endif
         hcla(i,jml-1-j)=hclanbt(i,j)*tex(i,jml-1-j,1)
	  enddo
      enddo
! avoid doublt counting
        jen = jml-jndg-4
      else
        jen = jml-2
      endif


! west
      if(ip_x .eq. 0) then
      do j=1,jml
      do i=1,indg+2
         hclawbt(i,j)=hclwbc(i,j,mb1)
     &        +(hclwbc(i,j,mb2)-hclwbc(i,j,mb1))*cbf
      enddo
      enddo

      do j=jsn,jen
      do i=3,indg+2
#ifdef BNDTIDE
       hclwbt(i,j) = h_tide(i+2,j)
#endif
       hcla(i+2,j)=(hcla(i+2,j)
     &        + chfhcl*c2dtts*(hclwbt(i,j)-hcl(i+2,j)))*tex(i+2,j,1)
       ! hcla(i+2,j)=hclwbt(i,j)*tex(i+2,j,1)
      enddo
      enddo

      do j=3,jml-2
      do i=1,2
#ifdef BNDTIDE
       hclawbt(i,j) = h_tide(i+2,j)
#endif
         hcla(i+2,j)=hclawbt(i,j)*tex(i+2,j,1)
      enddo
      enddo
      endif


!  east
      if(ip_x .eq. ipe-1) then
      do j=1,jml
      do i=1,indg+2
         hclaebt(i,j)=hclebc(i,j,mb1)
     $        +(hclebc(i,j,mb2)-hclebc(i,j,mb1))*cbf
      enddo
      enddo

      do j=jsn,jen
      do i=3,indg+2
#ifdef BNDTIDE
       hclebt(i,j) = h_tide(iml-1-i,j)
#endif
       hcla(iml-1-i,j)=(hcla(iml-1-i,j)
     &  + chfhcl*c2dtts*(hclebt(i,j)-hcl(iml-1-i,j)))*tex(iml-1-i,j,1)
        !hcla(iml-1-i,j)=hclebt(i,j)*tex(iml-1-i,j,1)
      enddo
      enddo

      do j=3,jml-2
      do i=1,2
#ifdef BNDTIDE
       hclaebt(i,j) = h_tide(iml-1-i,j)
#endif
         hcla(iml-1-i,j)=hclaebt(i,j)*tex(iml-1-i,j,1)
      enddo
      enddo
      endif
#endif
#ifdef BNDTIDE
      do i = 1,iml-1
      do j = 1,jml-1
         if(lw(ip_x)+i-1 >= 114 .and. lw(ip_x)+i-1 <= 160 .and.
     &         ls(ip_y)+j-1 >= 335) then
           hcla(i,j) = 0.
         endif
      enddo
      enddo
#endif


      do j = 2,jml-2
      do i = 2,iml-2
         hclu(i,j)=ex(i,j,1)*(ashf(j)*(hcl(i,j)+hcl(i+1,j))
     $        +anhf(j)*(hcl(i,j+1)+hcl(i+1,j+1)))/areauu(j)
         hclua(i,j)=ex(i,j,1)*(ashf(j)*(hcla(i,j)+hcla(i+1,j))
     $        +anhf(j)*(hcla(i,j+1)+hcla(i+1,j+1)))/areauu(j)
      enddo
      enddo

#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do j = 1,jml-1
         do i=1,2
          hclu(i+2,j)=(ashf(j)*(hclwbt(i,j)+hclwbt(i+1,j))
     &        +anhf(j)*(hclwbt(i,j+1)+hclwbt(i+1,j+1)))
     &        *ex(i+2,j,1)/areauu(j)
          hclua(i+2,j)=(ashf(j)*(hclawbt(i,j)+hclawbt(i+1,j))
     &        +anhf(j)*(hclawbt(i,j+1)+hclawbt(i+1,j+1)))
     &        *ex(i+2,j,1)/areauu(j)
         enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do j = 1,jml-1
         do i=1,2
          hclu(iml-2-i,j)=(ashf(j)*(hclebt(i,j)+hclebt(i+1,j))
     &        +anhf(j)*(hclebt(i,j+1)+hclebt(i+1,j+1)))
     &        *ex(iml-2-1,j,1)/areauu(j)
          hclua(iml-2-i,j)=(ashf(j)*(hclaebt(i,j)+hclaebt(i+1,j))
     &        +anhf(j)*(hclaebt(i,j+1)+hclaebt(i+1,j+1)))
     &        *ex(iml-2-1,j,1)/areauu(j)

         enddo
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do i=1,iml-1
         do j=1,2
          hclu(i,j+2)=(ashf(j+2)*(hclsbt(i,j)+hclsbt(i+1,j))
     &       +anhf(j+2)*(hclsbt(i,j+1)+hclsbt(i+1,j+1)))
     &        *ex(i,j+2,1)/areauu(j+2)
           hclua(i,j+2)=(ashf(j+2)*(hclasbt(i,j)+hclasbt(i+1,j))
     &        +anhf(j+2)*(hclasbt(i,j+1)+hclasbt(i+1,j+1)))
     &        *ex(i,j+2,1)/areauu(j+2)
         enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do i=1,iml-1
         do j=1,2
         hclu(i,jml-2-j)=(ashf(jml-2-j)*(hclnbt(i,j)+hclnbt(i+1,j))
     &        +anhf(jml-2-j)*(hclnbt(i,j+1)+hclnbt(i+1,j+1)))
     &        *ex(i,jml-2-j,1)/areauu(jml-2-j)
         hclua(i,jml-2-j)=(ashf(jml-2-j)*(hclanbt(i,j)+hclanbt(i+1,j))
     &        +anhf(jml-2-j)*(hclanbt(i,j+1)+hclanbt(i+1,j+1)))
     &        *ex(i,jml-2-j,1)/areauu(jml-2-j)
         enddo
      enddo
      endif
#endif

      call exch_t2d_n1p(hclu,hclua,1)
      call wait_t2d_n1p(hclu,hclua,1)
      call exch_t2d_e1p(hclu,hclua,55)
      call wait_t2d_e1p(hclu,hclua,55)

      do j = 2,jml-1
      do k = 1,kadp
      do i = 2,iml-1
         dzt(i,j,k) = tex(i,j,k)*(dz(k)+hcl(i,j)*dz(k)/depadp)
         volt(i,j,k) = tex(i,j,k)*(dz(k)
     $        + hcl(i,j)*dz(k)/depadp)*areat(i,j,k)
         voltr(i,j,k)=tex(i,j,k)/(volt(i,j,k)+(1.d0-tex(i,j,k)))

         dzu(i,j,k)=ex(i,j,k)*(dz(k)+hclu(i,j)*dz(k)/depadp)
         dzua(i,j,k)=ex(i,j,k)*(dz(k)+hclua(i,j)*dz(k)/depadp)
         volur(i,j,k)=ex(i,j,k)/(areauu(j)
     $        *(dz(k)+hclua(i,j)*dz(k)/depadp)+(1.d0-ex(i,j,k)))
      enddo
      enddo
      enddo

#ifdef BISMFRIC
      call exch_t3d_n2(dzu,35)
      call exch_t3d_s2(dzu,28)
      call wait_t3d_n2(dzu,35)
      call wait_t3d_s2(dzu,28)
      call exch_t3d_e2(dzu,70)
      call exch_t3d_w2(dzu,71)
#endif
c
      call wait_t3d_e1p(t,s,45,46)
      call wait_t3d_w1p(t,s,47,48)

      do k = 1,km
      do j = 2,jml-1
      do i = 2,iml-1
      sroot = dsqrt(s(i,j,k))
      rhoo(i,j,k) = 999.842594d0
     &  + t(i,j,k)*(6.793952d-2 + t(i,j,k)*(-9.095290d-3
     &    + t(i,j,k)*(1.001685d-4 + t(i,j,k)*(-1.120083d-6
     &    + 6.536332d-9*t(i,j,k)))))
     &  + s(i,j,k)*( (0.824493d0 + t(i,j,k)*(-4.0899d-3
     &    + t(i,j,k)*(7.6438d-5 + t(i,j,k)*(-8.2467d-7
     &    + 5.3875d-9*t(i,j,k))))) + sroot*(-5.72466d-3
     &    + t(i,j,k)*(1.0227d-4 - 1.6546d-6*t(i,j,k)))
     &    + 4.8314d-4*s(i,j,k))

      rhoo(i,j,k)=tex(i,j,k)*rhoo(i,j,k)*1.d-3-1.d0

!   compute rho(s,theta,p)

      stb0 = 1.965935d4
     &  + t(i,j,k)*(1.445863d2 + t(i,j,k)*(-1.722523d0
     &    + t(i,j,k)*(1.19238d-2 - 4.768276d-5*t(i,j,k))))
     &  + s(i,j,k)*((52.85624d0 + t(i,j,k)*(-3.1128126d-1
     &    + t(i,j,k)*(6.456036d-3 - 5.370396d-5*t(i,j,k))))
     &  + sroot*(3.884013d-1 + t(i,j,k)*(9.116446d-3
     &    - 4.628163d-4*t(i,j,k))))
      stb1 = 3.186518d0 + t(i,j,k)*(2.189412d-2
     &  + t(i,j,k)*(-2.823685d-4 + 1.715739d-6*t(i,j,k)))
     &  + s(i,j,k)*(6.703377d-3 + t(i,j,k)*(-1.839953d-4
     &    + t(i,j,k)*1.912264d-7) + 1.477291d-4*sroot)
      stb2 = 2.111102d-4 + t(i,j,k)*(-1.202016d-5
     &     + 1.364330d-7*t(i,j,k))
     &  + s(i,j,k)*(-2.048755d-6 + t(i,j,k)*(6.375979d-8*t(i,j,k)
     &    + 5.240967d-10*t(i,j,k)))

      stbm = stb0 + pd(k)*(stb1 + pd(k)*stb2)
      rho(i,j,k)=tex(i,j,k)
     &  *( (rhoo(i,j,k)+1.d0)/(1.d0 -pd(k)/stbm) -1.d0)

!   rhodh,rhouh is used in noh-mix
!   rhouf,rhodf,rhodh,rhouh is used in gm-dif

      stbdf = stb0 + pd(k+1)*(stb1 + pd(k+1)*stb2)
      rhodf(i,j,k) =tex(i,j,k)
     &   *( (rhoo(i,j,k)+1.d0)/(1.d0 -pd(k+1)/stbdf) -1.d0)

      stbdh = stb0 + pm(k+1)*(stb1 + pm(k+1)*stb2)
      rhodh(i,j,k)=tex(i,j,k)
     &   *( (rhoo(i,j,k)+1.d0)/(1.d0 -pm(k+1)/stbdh) -1.d0)

      stbuf = stb0 + pd(k-1)*(stb1 + pd(k-1)*stb2)
      rhouf(i,j,k)=tex(i,j,k)
     &   *( (rhoo(i,j,k)+1.d0)/(1.d0 -pd(k-1)/stbuf) -1.d0)

      stbuh = stb0 + pm(k)*(stb1 + pm(k)*stb2)
      rhouh(i,j,k)=tex(i,j,k)
     &   *( (rhoo(i,j,k)+1.d0)/(1.d0 -pm(k)/stbuh) -1.d0)

      enddo
      enddo
      enddo

      do j = 3,jml-2
      do i = 3,iml-2
ccc   bottom stress
ccc   k is the bottom level, but in land k=1 for use in array subscript

      k=nint(ex(i,j,1))*exn(i,j)+1-nint(ex(i,j,1))
      botvel=ex(i,j,k)*areauu(j)*1.225d-3 ! Ishikawa (original model)
      ! botvel=ex(i,j,k)*areauu(j)*2.6d-3   ! Nakamura et al. 2000, Hatayama et al. 1996 (barotropic model)
      ! botvel=ex(i,j,k)*areauu(j)*1.6d-3   ! Hu et al. 2010 (0.0015-0.0017) (barotropic model)
      ! botvel=ex(i,j,k)*areauu(j)*2.0d-3     ! Schiller 2004, Schiller and Fiedler 2007 (MOM4 model)
     $     *dsqrt(ub(i,j,k)*ub(i,j,k)+vb(i,j,k)*vb(i,j,k))
      ua(i,j,k)=ua(i,j,k)-botvel*(ub(i,j,k)*bcs-vb(i,j,k)*sgn(j))
      va(i,j,k)=va(i,j,k)-botvel*(vb(i,j,k)*bcs+ub(i,j,k)*sgn(j))

!cc   wind stress
#ifdef BNDTIDE
       wsx(i,j) = 0.
       wsy(i,j) = 0.
       netq_o(i,j) = 0.
       wflux(i,j) = 0.
#endif
      ua(i,j,1)=ua(i,j,1)+ex(i,j,1)*wsx(i,j)*factws
     &     *areauu(j)*ddmnar
      va(i,j,1)=va(i,j,1)+ex(i,j,1)*wsy(i,j)*factws
     $     *areauu(j)*ddmnar

#ifdef EQ100M
      ua(i,j,1)=ua(i,j,1)+gref2(i,j,1)*(uref-u(i,j,1))
     &     *areauu(j)*(dz(k)+hclua(i,j)*dz(k)/depadp)
      va(i,j,1)=va(i,j,1)+gref2(i,j,1)*(vref-v(i,j,1))
     &     *areauu(j)*(dz(k)+hclua(i,j)*dz(k)/depadp)
#endif
c
ccc   temp surface flux including water flux
ccc   rain temperature is same as SST now
c
      ta(i,j,1)=ta(i,j,1)+tex(i,j,1)*areat(i,j,1)*(
     &     netq_o(i,j)*0.24d0*ddmnar+t(i,j,1)*wflux(i,j))
      enddo
      enddo

!   surface flux and body forcing
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
         gref1=tex(i,j,k)*gref(i,j,k)
#ifdef PC68M
         gref1=tex(i,j,k)
     &        *dmax1( dble(gref(i,j,k)),dble(grefini(j))*fcgini )
#endif
#ifdef NWNPAC
#if defined(CASE58) || defined(CASE60)
         gref1=tex(i,j,k)
     &        *dmax1( dble(gref(i,j,k)),dble(grefini(j))*fcgini )
#else
         gref1=tex(i,j,k)*gref(i,j,k)*fcgini
#endif
#endif
         ta(i,j,k) = ta(i,j,k) +
     $        tex(i,j,k)*gref1*(tref(i,j,k)-t(i,j,k))
     &        *areat(i,j,k)*dzt(i,j,k)
         sa(i,j,k) = sa(i,j,k) +
     $        tex(i,j,k)*gref1*(sref(i,j,k)-s(i,j,k))
     &        *areat(i,j,k)*dzt(i,j,k)
      enddo
      enddo
      enddo

!   presure gradient
      do j = 3,jml-2

!   k = 1
      do i = 3,iml-2
         dpdx(i,j,1)=dpdx(i,j,0)
     &        +ex(i,j,1)*grav*dx2r*csr(j)*ddmnar*0.5d0*dzu(i,j,1)
     $        *(rho(i+1,j+1,1)+rho(i+1,j,1)-rho(i,j+1,1)-rho(i,j,1))
         dpdy(i,j,1)=dpdy(i,j,0)
     &        +ex(i,j,1)*grav*dy2r*ddmnar*0.5d0*dzu(i,j,1)
     $        *(rho(i+1,j+1,1)-rho(i+1,j,1)+rho(i,j+1,1)-rho(i,j,1))
      enddo

      do k = 2,kadp
      do i = 3,iml-2
         dpdx(i,j,k)=dpdx(i,j,k-1)+ex(i,j,k)
     $        *grav*ddmnar*0.5d0*(dzu(i,j,k)+dzu(i,j,k-1))
     $        *dx2r*( 0.5*csr(j)*(rho(i+1,j+1,k-1)+rho(i+1,j,k-1)
     $        -rho(i,j+1,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i+1,j,k)-rho(i,j+1,k)-rho(i,j,k))
     $        +(hcl(i+1,j+1)-hcl(i,j+1))/depadp*0.25d0
     $        *cstr(j+1)*(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j+1,k)-rho(i,j+1,k))
     $        +(hcl(i+1,j)-hcl(i,j))/depadp*0.25d0
     $        *cstr(j)*(rho(i+1,j,k-1)+rho(i,j,k-1)
     $        -rho(i+1,j,k)-rho(i,j,k)) )
         dpdy(i,j,k)=dpdy(i,j,k-1)+ex(i,j,k)
     $        *grav*ddmnar*.5d0*(dzu(i,j,k)+dzu(i,j,k-1))
     $        *dy2r*( 0.5d0*(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i,j+1,k)-rho(i+1,j,k)-rho(i,j,k))
     $        +(hcl(i+1,j+1)+hcl(i,j+1)-hcl(i+1,j)-hcl(i,j))
     &        /depadp*0.25d0
     $        *(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)+rho(i+1,j,k-1)
     $        +rho(i,j,k-1)-rho(i+1,j+1,k)-rho(i,j+1,k)
     $        -rho(i+1,j,k)-rho(i,j,k)) )
      enddo
      enddo

      k = kadp+1

      do i = 3,iml-2
         dpdx(i,j,k)=dpdx(i,j,k-1)+ex(i,j,k)*csr(j)
     $        *grav*dx4r*ddmnar
     $        *(dzz(k)*.5d0*(1.d0+dzr(k)*dzu(i,j,k))
     $        +0.5d0*hclu(i,j)*dep(kadp)/depadp)
     $        *(rho(i+1,j+1,k-1)+rho(i+1,j,k-1)
     $        -rho(i,j+1,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i+1,j,k)-rho(i,j+1,k)-rho(i,j,k))
         dpdy(i,j,k)=dpdy(i,j,k-1)+ex(i,j,k)
     $        *grav*dy4r*ddmnar
     $        *(dzz(k)*.5d0*(1.d0+dzr(k)*dzu(i,j,k))
     $        +0.5d0*hclu(i,j)*dep(kadp)/depadp)
     $        *(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i,j+1,k)-rho(i+1,j,k)-rho(i,j,k))
      enddo

      do k = kadp+2,km
      do i = 3,iml-2
         dpdx(i,j,k)=dpdx(i,j,k-1)+ex(i,j,k)*csr(j)
     $        *grav*dx4r*dzz(k)*ddmnar*.5d0*(1.d0+dzr(k)*dzu(i,j,k))
     $        *(rho(i+1,j+1,k-1)+rho(i+1,j,k-1)
     $        -rho(i,j+1,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i+1,j,k)-rho(i,j+1,k)-rho(i,j,k))
         dpdy(i,j,k)=dpdy(i,j,k-1)+ex(i,j,k)
     $        *grav*dy4r*dzz(k)*ddmnar*.5d0*(1.d0+dzr(k)*dzu(i,j,k))
     $        *(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i,j+1,k)-rho(i+1,j,k)-rho(i,j,k))
      enddo
      enddo

      enddo

#ifdef BNDTIDE
       dpdx(:,:,:) = 0.
       dpdy(:,:,:) = 0.
#endif

!   definition of ustar and vstar

      call wait_t3d_e1p(u,v,78,79)
      call wait_t3d_w1p(u,v,56,57)

      do k = 1,km
      do j = 3,jml-1
      do i = 3,iml-1
         ustar(i,j,k)=0.5d0*dy*((2.d0-ex(i,j-1,k))*u(i,j,k)*dzu(i,j,k)
     &        +(2.d0-ex(i,j,k))*u(i,j-1,k)*dzu(i,j-1,k))
         vstar(i,j,k)=0.5d0*dx*((2.d0-ex(i-1,j,k))*v(i,j,k)*dzu(i,j,k)
     &        +(2.d0-ex(i,j,k))*v(i-1,j,k)*dzu(i-1,j,k))*cs(j)
      enddo
      enddo
      enddo

      call exch_t3d_s1(vstar,24)
      call exch_t3d_w1(ustar,80)

      call wait_t3d_e2p(tb,sb,49,50)
      call wait_t3d_e1p(tb,sb,51,52)
      call wait_t3d_w1p(tb,sb,53,54)

!   horizontal advection and horizontal diffusion of t and s

      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
c
ccc   upwind and centered
c
c        if(ustar(i,j,k).gt.0.d0) then
c          facw = (1.d0+hupp)*0.5d0
c          face = (1.d0-hupp)*0.5d0
c        else
c          facw = (1.d0-hupp)*0.5d0
c          face = (1.d0+hupp)*0.5d0
c        endif
cc
c        if(vstar(i,j,k).gt.0.d0) then
c          facs = (1.d0+hupp)*0.5d0
c          facn = (1.d0-hupp)*0.5d0
c        else
c          facs = (1.d0-hupp)*0.5d0
c          facn = (1.d0+hupp)*0.5d0
c        endif
cc
c        tust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
c     $        *(facw*t(i,j,k)+face*t(i+1,j,k))
c        sust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
c     $        *(facw*s(i,j,k)+face*s(i+1,j,k))
c
c        tvst(i,j,k)=.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
c     $       *(facs*t(i,j,k)+facn*t(i,j+1,k))
c        svst(i,j,k)=.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
c     $       *(facs*s(i,j,k)+facn*s(i,j+1,k))
c
ccc   in case of utopia
c
c         tee = tex(i+2,j,k)*tb(i+2,j,k)
c     &        +(1.-tex(i+2,j,k))*tb(i+1,j,k)
c         tww = tex(i-1,j,k)*tb(i-1,j,k)
c     &        +(1.-tex(i-1,j,k))*tb(i,j,k)
c         tnn = tex(i,j+2,k)*tb(i,j+2,k)
c     &        +(1.-tex(i,j+2,k))*tb(i,j+1,k)
c         tss = tex(i,j-1,k)*tb(i,j-1,k)
c     &        +(1.-tex(i,j-1,k))*tb(i,j,k)
c         see = tex(i+2,j,k)*sb(i+2,j,k)
c     &        +(1.-tex(i+2,j,k))*sb(i+1,j,k)
c         sww = tex(i-1,j,k)*sb(i-1,j,k)
c     &        +(1.-tex(i-1,j,k))*sb(i,j,k)
c         snn = tex(i,j+2,k)*sb(i,j+2,k)
c     &        +(1.-tex(i,j+2,k))*sb(i,j+1,k)
c         sss = tex(i,j-1,k)*sb(i,j-1,k)
c     &        +(1.-tex(i,j-1,k))*sb(i,j,k)
c
         texex=tex(i+2,j,k)*(ex(i+1,j-1,k)+ex(i+1,j,k)
     &        -ex(i+1,j-1,k)*ex(i+1,j,k))
         tee=texex*tb(i+2,j,k)+(1.-texex)*tb(i+1,j,k)
         see=texex*sb(i+2,j,k)+(1.-texex)*sb(i+1,j,k)
         texex=tex(i-1,j,k)*(ex(i-1,j-1,k)+ex(i-1,j,k)
     &        -ex(i-1,j-1,k)*ex(i-1,j,k))
         tww=texex*tb(i-1,j,k)+(1.-texex)*tb(i,j,k)
         sww=texex*sb(i-1,j,k)+(1.-texex)*sb(i,j,k)
         texex=tex(i,j+2,k)*(ex(i-1,j+1,k)+ex(i,j+1,k)
     &        -ex(i-1,j+1,k)*ex(i,j+1,k))
         tnn=texex*tb(i,j+2,k)+(1.-texex)*tb(i,j+1,k)
         snn=texex*sb(i,j+2,k)+(1.-texex)*sb(i,j+1,k)
         texex=tex(i,j-1,k)*(ex(i-1,j-1,k)+ex(i,j-1,k)
     &        -ex(i-1,j-1,k)*ex(i,j-1,k))
         tss=texex*tb(i,j-1,k)+(1.-texex)*tb(i,j,k)
         sss=texex*sb(i,j-1,k)+(1.-texex)*sb(i,j,k)

        cue = 0.5d0*((2.d0-ex(i,j-1,k))*u(i,j,k)
     &    +(2.d0-ex(i,j,k))*u(i,j-1,k))*c2dtts
        cve = 0.5d0*((2.d0-ex(i,j-1,k))*v(i,j,k)
     &    +(2.d0-ex(i,j,k))*v(i,j-1,k))*c2dtts
        cun=0.5d0*((2.d0-ex(i-1,j,k))*u(i,j,k)
     &    +(2.d0-ex(i,j,k))*u(i-1,j,k) )*c2dtts
        cvn=0.5d0*((2.d0-ex(i-1,j,k))*v(i,j,k)
     &    +(2.d0-ex(i,j,k))*v(i-1,j,k) )*c2dtts

        if(ustar(i,j,k).gt.0.d0) then
          tc20 = 0.5*(tb(i+1,j,k)-2.*tb(i,j,k)+tww)
     &           *dxr*dxr*cstr(j)*cstr(j)
          tc02 = 0.5*( (tb(i,j+1,k)-tb(i,j,k))*tex(i,j+1,k)
     &         *(ex(i-1,j,k)+ex(i,j,k)-ex(i-1,j,k)*ex(i,j,k))
     &         -(tb(i,j,k)-tb(i,j-1,k))*tex(i,j-1,k)
     &         *(ex(i-1,j-1,k)+ex(i,j-1,k)-ex(i-1,j-1,k)*ex(i,j-1,k))
     &         )*dyr*dyr
          sc20 = 0.5*(sb(i+1,j,k)-2.*sb(i,j,k)+sww)
     &           *dxr*dxr*cstr(j)*cstr(j)
          sc02 = 0.5*( (sb(i,j+1,k)-sb(i,j,k))*tex(i,j+1,k)
     &         *(ex(i-1,j,k)+ex(i,j,k)-ex(i-1,j,k)*ex(i,j,k))
     &         -(sb(i,j,k)-sb(i,j-1,k))*tex(i,j-1,k)
     &         *(ex(i-1,j-1,k)+ex(i,j-1,k)-ex(i-1,j-1,k)*ex(i,j-1,k))
     &         )*dyr*dyr
        else
          tc20 = 0.5*(tee - 2.*tb(i+1,j,k)+tb(i,j,k))
     &           *dxr*dxr*cstr(j)*cstr(j)
          tc02 = 0.5*( (tb(i+1,j+1,k)-tb(i+1,j,k))*tex(i+1,j+1,k)
     &         *(ex(i,j,k)+ex(i+1,j,k)-ex(i,j,k)*ex(i+1,j,k))
     &         -(tb(i+1,j,k)-tb(i+1,j-1,k))*tex(i+1,j-1,k)
     &         *(ex(i,j-1,k)+ex(i+1,j-1,k)-ex(i,j-1,k)*ex(i+1,j-1,k))
     &         )*dyr*dyr
          sc20 = 0.5*(see - 2.*sb(i+1,j,k)+sb(i,j,k))
     &           *dxr*dxr*cstr(j)*cstr(j)
          sc02 = 0.5*( (sb(i+1,j+1,k)-sb(i+1,j,k))*tex(i+1,j+1,k)
     &         *(ex(i,j,k)+ex(i+1,j,k)-ex(i,j,k)*ex(i+1,j,k))
     &         -(sb(i+1,j,k)-sb(i+1,j-1,k))*tex(i+1,j-1,k)
     &         *(ex(i,j-1,k)+ex(i+1,j-1,k)-ex(i,j-1,k)*ex(i+1,j-1,k))
     &         )*dyr*dyr
        endif

        if(cve.gt.0.d0) then
          tc11 = dyr*dxr*cstr(j)*
     &          ((tb(i,j-1,k)-tb(i,j,k))*tex(i,j-1,k)
     &          *(ex(i-1,j-1,k)+ex(i,j-1,k)-ex(i-1,j-1,k)*ex(i,j-1,k))
     &          -(tb(i+1,j-1,k)-tb(i+1,j,k))*tex(i+1,j-1,k)
     &          *(ex(i,j-1,k)+ex(i+1,j-1,k)-ex(i,j-1,k)*ex(i+1,j-1,k)))
          sc11 = dyr*dxr*cstr(j)*
     &         ((sb(i,j-1,k)-sb(i,j,k))*tex(i,j-1,k)
     &         *(ex(i-1,j-1,k)+ex(i,j-1,k)-ex(i-1,j-1,k)*ex(i,j-1,k))
     &         -(sb(i+1,j-1,k)-sb(i+1,j,k))*tex(i+1,j-1,k)
     &         *(ex(i,j-1,k)+ex(i+1,j-1,k)-ex(i,j-1,k)*ex(i+1,j-1,k)))
        else
          tc11 = dxr*dyr*cstr(j+1)
     &          *((tb(i+1,j+1,k)-tb(i+1,j,k))*tex(i+1,j+1,k)
     &          *(ex(i,j,k)+ex(i+1,j,k)-ex(i,j,k)*ex(i+1,j,k))
     &          -(tb(i,j+1,k)-tb(i,j,k))*tex(i,j+1,k)
     &          *(ex(i-1,j,k)+ex(i,j,k)-ex(i-1,j,k)*ex(i,j,k)))
          sc11 = dxr*dyr*cstr(j+1)
     &         *((sb(i+1,j+1,k)-sb(i+1,j,k))*tex(i+1,j+1,k)
     &         *(ex(i,j,k)+ex(i+1,j,k)-ex(i,j,k)*ex(i+1,j,k))
     &         -(sb(i,j+1,k)-sb(i,j,k))*tex(i,j+1,k)
     &         *(ex(i-1,j,k)+ex(i,j,k)-ex(i-1,j,k)*ex(i,j,k)))
        endif

        if(ustar(i,j,k)*cve.gt.0.d0) then
           tc01 = ( (tb(i,j+1,k)-tb(i,j,k))*tex(i,j+1,k)
     &          *(ex(i-1,j,k)+ex(i,j,k)-ex(i-1,j,k)*ex(i,j,k))
     &          +tex(i+1,j-1,k)*(tb(i+1,j,k)-tb(i+1,j-1,k))
     &          *(ex(i,j-1,k)+ex(i+1,j-1,k)-ex(i,j-1,k)*ex(i+1,j-1,k))
     &          )*dyr
           sc01 = ( (sb(i,j+1,k)-sb(i,j,k))*tex(i,j+1,k)
     &          *(ex(i-1,j,k)+ex(i,j,k)-ex(i-1,j,k)*ex(i,j,k))
     &          +tex(i+1,j-1,k)*(sb(i+1,j,k)-sb(i+1,j-1,k))
     &          *(ex(i,j-1,k)+ex(i+1,j-1,k)-ex(i,j-1,k)*ex(i+1,j-1,k))
     &          )*dyr
        else
           tc01 = ( (tb(i+1,j+1,k)-tb(i+1,j,k))*tex(i+1,j+1,k)
     &          *(ex(i,j,k)+ex(i+1,j,k)-ex(i,j,k)*ex(i+1,j,k))
     &          +tex(i,j-1,k)*(tb(i,j,k)-tb(i,j-1,k))
     &          *(ex(i-1,j-1,k)+ex(i,j-1,k)-ex(i-1,j-1,k)*ex(i,j-1,k))
     &          )*dyr
           sc01 = ( (sb(i+1,j+1,k)-sb(i+1,j,k))*tex(i+1,j+1,k)
     &          *(ex(i,j,k)+ex(i+1,j,k)-ex(i,j,k)*ex(i+1,j,k))
     &          +tex(i,j-1,k)*(sb(i,j,k)-sb(i,j-1,k))
     &          *(ex(i-1,j-1,k)+ex(i,j-1,k)-ex(i-1,j-1,k)*ex(i,j-1,k))
     &          )*dyr
        endif

        if(vstar(i,j,k).gt.0.d0) then
           td02 = 0.5*(tb(i,j+1,k)-2.*tb(i,j,k)+tss)*dyr*dyr
           td20 = 0.5*( (tb(i+1,j,k)-tb(i,j,k))*tex(i+1,j,k)
     &          *(ex(i,j-1,k)+ex(i,j,k)-ex(i,j-1,k)*ex(i,j,k))
     &          -(tb(i,j,k)-tb(i-1,j,k))*tex(i-1,j,k)
     &          *(ex(i-1,j-1,k)+ex(i-1,j,k)-ex(i-1,j-1,k)*ex(i-1,j,k))
     &          )*dxr*dxr*cstr(j)*cstr(j)
           sd02 = 0.5*(sb(i,j+1,k)-2.*sb(i,j,k)+sss)*dyr*dyr
           sd20 = 0.5*((sb(i+1,j,k)-sb(i,j,k))*tex(i+1,j,k)
     &          *(ex(i,j-1,k)+ex(i,j,k)-ex(i,j-1,k)*ex(i,j,k))
     &          -(sb(i,j,k)-sb(i-1,j,k))*tex(i-1,j,k)
     &          *(ex(i-1,j-1,k)+ex(i-1,j,k)-ex(i-1,j-1,k)*ex(i-1,j,k))
     &          )*dxr*dxr*cstr(j)*cstr(j)
        else
           td02 = 0.5*(tnn-2.*tb(i,j+1,k)+tb(i,j,k))*dyr*dyr
           td20 = 0.5*((tb(i+1,j+1,k)-tb(i,j+1,k))*tex(i+1,j+1,k)
     &          *(ex(i,j,k)+ex(i,j+1,k)-ex(i,j,k)*ex(i,j+1,k))
     &          -(tb(i,j+1,k)-tb(i-1,j+1,k))*tex(i-1,j+1,k)
     &          *(ex(i-1,j,k)+ex(i-1,j+1,k)-ex(i-1,j,k)*ex(i-1,j+1,k))
     &          )*dxr*dxr*cstr(j+1)*cstr(j+1)
           sd02 = 0.5*(snn-2.*sb(i,j+1,k)+sb(i,j,k))*dyr*dyr
           sd20 = 0.5*((sb(i+1,j+1,k)-sb(i,j+1,k))*tex(i+1,j+1,k)
     &          *(ex(i,j,k)+ex(i,j+1,k)-ex(i,j,k)*ex(i,j+1,k))
     &          -(sb(i,j+1,k)-sb(i-1,j+1,k))*tex(i-1,j+1,k)
     &          *(ex(i-1,j,k)+ex(i-1,j+1,k)-ex(i-1,j,k)*ex(i-1,j+1,k))
     &          )*dxr*dxr*cstr(j+1)*cstr(j+1)
        endif

        if(cun.gt.0.d0) then
           td11 = 0.5*((tb(i-1,j,k)-tb(i,j,k))*tex(i-1,j,k)
     &          *(ex(i-1,j-1,k)+ex(i-1,j,k)-ex(i-1,j-1,k)*ex(i-1,j,k))
     &          -(tb(i-1,j+1,k)-tb(i,j+1,k))*tex(i-1,j+1,k)
     &          *(ex(i-1,j,k)+ex(i-1,j+1,k)-ex(i-1,j,k)*ex(i-1,j+1,k))
     &          )*dxr*cstr(j)*dyr
           sd11 = 0.5*((sb(i-1,j,k)-sb(i,j,k))*tex(i-1,j,k)
     &          *(ex(i-1,j-1,k)+ex(i-1,j,k)-ex(i-1,j-1,k)*ex(i-1,j,k))
     &          -(sb(i-1,j+1,k)-sb(i,j+1,k))*tex(i-1,j+1,k)
     &          *(ex(i-1,j,k)+ex(i-1,j+1,k)-ex(i-1,j,k)*ex(i-1,j+1,k))
     &          )*dxr*cstr(j)*dyr
        else
           td11 = 0.5*((tb(i+1,j+1,k)-tb(i,j+1,k))*tex(i+1,j+1,k)
     &          *(ex(i,j,k)+ex(i,j+1,k)-ex(i,j,k)*ex(i,j+1,k))
     &          -(tb(i+1,j,k)-tb(i,j,k))*tex(i+1,j,k)
     &          *(ex(i,j-1,k)+ex(i,j,k)-ex(i,j-1,k)*ex(i,j,k))
     &          )*dxr*cstr(j)*dyr
           sd11 = 0.5*((sb(i+1,j+1,k)-sb(i,j+1,k))*tex(i+1,j+1,k)
     &          *(ex(i,j,k)+ex(i,j+1,k)-ex(i,j,k)*ex(i,j+1,k))
     &          -(sb(i+1,j,k)-sb(i,j,k))*tex(i+1,j,k)
     &          *(ex(i,j-1,k)+ex(i,j,k)-ex(i,j-1,k)*ex(i,j,k))
     &          )*dxr*cstr(j)*dyr
        endif

        if(cun*vstar(i,j,k).gt.0.d0) then
           td10 = 0.5*((tb(i+1,j,k)-tb(i,j,k))*tex(i+1,j,k)
     &          *(ex(i,j-1,k)+ex(i,j,k)-ex(i,j-1,k)*ex(i,j,k))
     &          -(tb(i-1,j+1,k)-tb(i,j+1,k))*tex(i-1,j+1,k)
     &          *(ex(i,j,k)+ex(i,j+1,k)-ex(i,j,k)*ex(i,j+1,k))
     &          )*dxr*cstr(j)
           sd10 = 0.5*((sb(i+1,j,k)-sb(i,j,k))*tex(i+1,j,k)
     &          *(ex(i,j-1,k)+ex(i,j,k)-ex(i,j-1,k)*ex(i,j,k))
     &       -(sb(i-1,j+1,k)-sb(i,j+1,k))*tex(i-1,j+1,k)
     &          *(ex(i,j,k)+ex(i,j+1,k)-ex(i,j,k)*ex(i,j+1,k))
     &           )*dxr*cstr(j)
        else
           td10 = 0.5*((tb(i+1,j+1,k)-tb(i,j+1,k))*tex(i+1,j+1,k)
     &          *(ex(i+1,j,k)+ex(i+1,j+1,k)-ex(i+1,j,k)*ex(i+1,j+1,k))
     &          -(tb(i-1,j,k)-tb(i,j,k))*tex(i-1,j,k)
     &          *(ex(i-1,j-1,k)+ex(i-1,j,k)-ex(i-1,j-1,k)*ex(i-1,j,k))
     &          )*dxr*cstr(j)
           sd10 = 0.5*((sb(i+1,j+1,k)-sb(i,j+1,k))*tex(i+1,j+1,k)
     &          *(ex(i+1,j,k)+ex(i+1,j+1,k)-ex(i+1,j,k)*ex(i+1,j+1,k))
     &          -(sb(i-1,j,k)-sb(i,j,k))*tex(i-1,j,k)
     &          *(ex(i-1,j-1,k)+ex(i-1,j,k)-ex(i-1,j-1,k)*ex(i-1,j,k))
     &          )*dxr*cstr(j)
        endif

        tc00 = 0.5*(tb(i,j,k)+tb(i+1,j,k))
     &    -0.25*tc20*dx*dx*cst(j)*cst(j)
        tc10 = (tb(i+1,j,k)-tb(i,j,k))*dxr*cstr(j)
        td00 = 0.5*(tb(i,j,k)+tb(i,j+1,k))
     &   -0.25*td02*dy*dy
        td01 = (tb(i,j+1,k)-tb(i,j,k))*dyr
        sc00 = 0.5*(sb(i,j,k)+sb(i+1,j,k))
     &    -0.25*sc20*dx*dx*cst(j)*cst(j)
        sc10 = (sb(i+1,j,k)-sb(i,j,k))*dxr*cstr(j)
        sd00 = 0.5*(sb(i,j,k)+sb(i,j+1,k))
     &    -0.25*sd02*dy*dy
        sd01 = (sb(i,j+1,k)-sb(i,j,k))*dyr

        tust(i,j,k) = 0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &      *(tc00 - 0.5*cue*tc10
     &    +(cue*cue-0.25*dx*dx*cst(j)*cst(j))/3.*tc20
     &       -0.5*cve*tc01+0.5*cve*cve*tc02 + 0.25*cue*cve*tc11)
        sust(i,j,k) = 0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &     *(sc00 - 0.5*cue*sc10
     &    +(cue*cue-0.25*dx*dx*cst(j)*cst(j))/3.*sc20
     &       -0.5*cve*sc01+0.5*cve*cve*sc02 + 0.25*cue*cve*sc11)

        tvst(i,j,k) = .5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &    *(td00 - 0.5*cvn*td01
     &    +(cve*cve-0.25*dy*dy)/3.*td02 - 0.5*cun*td10
     &    +0.5*cun*cun*td20 + 0.25*cun*cvn*td11)
        svst(i,j,k) = .5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &    *(sd00 - 0.5*cvn*sd01
     &    +(cve*cve-0.25*dy*dy)/3.*sd02 - 0.5*cun*sd10
     &    +0.5*cun*cun*sd20 + 0.25*cun*cvn*sd11)

!   off-shore flow is downstream
        if(ustar(i,j,k).gt.0..and.tex(i-1,j,k).eq.0.)then
           tust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *tb(i,j,k)
           sust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *sb(i,j,k)
        endif
        if(ustar(i,j,k).lt.0..and.tex(i+2,j,k).eq.0.)then
           tust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *tb(i+1,j,k)
           sust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *sb(i+1,j,k)
        endif
        if(vstar(i,j,k).gt.0..and.tex(i,j-1,k).eq.0.)then
           tvst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *tb(i,j,k)
           svst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *sb(i,j,k)
        endif
        if(vstar(i,j,k).lt.0..and.tex(i,j+2,k).eq.0.)then
           tvst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *tb(i,j+1,k)
           svst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *sb(i,j+1,k)
        endif
!   bottom flow is downstream
        if(ex(i-1,j-1,k+1)*ex(i,j-1,k+1)*ex(i-1,j,k+1)*ex(i,j,k+1)
     &       .eq.0.)then
        if(ustar(i,j,k).gt.0.)then
           tust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *tb(i,j,k)
           sust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *sb(i,j,k)
        endif
        if(ustar(i,j,k).lt.0.)then
           tust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *tb(i+1,j,k)
           sust(i,j,k)=0.5d0*ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &          *sb(i+1,j,k)
        endif
        if(vstar(i,j,k).gt.0.)then
           tvst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *tb(i,j,k)
           svst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *sb(i,j,k)
        endif
        if(vstar(i,j,k).lt.0.)then
           tvst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *tb(i,j+1,k)
           svst(i,j,k)=0.5d0*vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &          *sb(i,j+1,k)
        endif
        endif

ccc   normal diffusion > see gm-dif
c
c        tudf(i,j,k)=
c     $   -0.5d0*dydxr*hdts*cstr(j)*(dzu(i,j,k)+dzu(i,j-1,k))
c     $         *(tb(i+1,j,k)-tb(i,j,k))
c     &       *tex(i,j,k)*tex(i+1,j,k)
c        sudf(i,j,k)=
c     $   -0.5d0*dydxr*hdts*cstr(j)*(dzu(i,j,k)+dzu(i,j-1,k))
c     $         *(sb(i+1,j,k)-sb(i,j,k))
c     &       *tex(i,j,k)*tex(i+1,j,k)
c        tvdf(i,j,k)=
c     $    -.5d0*dxdyr*hdts*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
c     $       *(tb(i,j+1,k)-tb(i,j,k))
c     &       *tex(i,j,k)*tex(i,j+1,k)
c        svdf(i,j,k)=
c     $    -.5d0*dxdyr*hdts*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
c     $       *(sb(i,j+1,k)-sb(i,j,k))
c     &       *tex(i,j,k)*tex(i,j+1,k)
c
      enddo
      enddo
      enddo

      call exch_t3d_s1p(tvst,svst,4,5)
      call exch_t3d_w1p(tust,sust,76,77)

      call wait_t3d_s1(vstar,24)
      call wait_t3d_w1(ustar,80)

      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
         uadv=(
     $    (u(i,j,k)+u(i-1,j,k))*(
     $      cxn(i,j,k)*(ustar(i-1,j,k)+ustar(i,j,k))
     $     +cxs(i,j+1,k)*(ustar(i-1,j+1,k)+ustar(i,j+1,k)) )
     &   -(u(i,j,k)+u(i+1,j,k))*(
     $      cxn(i+1,j,k)*(ustar(i,j,k)+ustar(i+1,j,k))
     $     +cxs(i+1,j+1,k)*(ustar(i,j+1,k)+ustar(i+1,j+1,k)) )
     $   +(u(i,j,k)+u(i,j-1,k))*(
     $      cye(i,j,k)*(vstar(i,j-1,k)+vstar(i,j,k))
     $     +cyw(i+1,j,k)*(vstar(i+1,j-1,k)+vstar(i+1,j,k)) )
     $   -(u(i,j,k)+u(i,j+1,k))*(
     $      cye(i,j+1,k)*(vstar(i,j,k)+vstar(i,j+1,k))
     $     +cyw(i+1,j+1,k)*(vstar(i+1,j,k)+vstar(i+1,j+1,k)) )
     $   +(u(i,j,k)+u(i-1,j-1,k))*(
     $      cne(i,j,k)*(ustar(i-1,j,k)+ustar(i,j,k)
     $       +vstar(i,j-1,k)+vstar(i,j,k)) )
     $   -(u(i,j,k)+u(i+1,j+1,k))*(
     $      cne(i+1,j+1,k)*(ustar(i,j+1,k)+ustar(i+1,j+1,k)
     $       +vstar(i+1,j,k)+vstar(i+1,j+1,k)) )
     $   +(u(i-1,j+1,k)+u(i,j,k))*(
     $      cse(i,j+1,k)*(ustar(i-1,j+1,k)+ustar(i,j+1,k)
     $       -vstar(i,j,k)-vstar(i,j+1,k)) )
     $   -(u(i,j,k)+u(i+1,j-1,k))*(
     $      cse(i+1,j,k)*(ustar(i,j,k)+ustar(i+1,j,k)
     $       -vstar(i+1,j-1,k)-vstar(i+1,j,k)) ))/24.
     $       +u(i,j,k)*v(i,j,k)*tng(j)/radius*dzu(i,j,k)*areauu(j)
        vadv=(
     $    (v(i,j,k)+v(i-1,j,k))*(
     $      cxn(i,j,k)*(ustar(i-1,j,k)+ustar(i,j,k))
     $     +cxs(i,j+1,k)*(ustar(i-1,j+1,k)+ustar(i,j+1,k)) )
     $   -(v(i,j,k)+v(i+1,j,k))*(
     $      cxn(i+1,j,k)*(ustar(i,j,k)+ustar(i+1,j,k))
     $     +cxs(i+1,j+1,k)*(ustar(i,j+1,k)+ustar(i+1,j+1,k)) )
     $   +(v(i,j,k)+v(i,j-1,k))*(
     $      cye(i,j,k)*(vstar(i,j-1,k)+vstar(i,j,k))
     $     +cyw(i+1,j,k)*(vstar(i+1,j-1,k)+vstar(i+1,j,k)) )
     $   -(v(i,j,k)+v(i,j+1,k))*(
     $      cye(i,j+1,k)*(vstar(i,j,k)+vstar(i,j+1,k))
     $     +cyw(i+1,j+1,k)*(vstar(i+1,j,k)+vstar(i+1,j+1,k)) )
     $   +(v(i,j,k)+v(i-1,j-1,k))*(
     $      cne(i,j,k)*(ustar(i-1,j,k)+ustar(i,j,k)
     $       +vstar(i,j-1,k)+vstar(i,j,k)) )
     $   -(v(i,j,k)+v(i+1,j+1,k))*(
     $      cne(i+1,j+1,k)*(ustar(i,j+1,k)+ustar(i+1,j+1,k)
     $       +vstar(i+1,j,k)+vstar(i+1,j+1,k)) )
     $   +(v(i-1,j+1,k)+v(i,j,k))*(
     $      cse(i,j+1,k)*(ustar(i-1,j+1,k)+ustar(i,j+1,k)
     $       -vstar(i,j,k)-vstar(i,j+1,k)) )
     $   -(v(i,j,k)+v(i+1,j-1,k))*(
     $      cse(i+1,j,k)*(ustar(i,j,k)+ustar(i+1,j,k)
     $       -vstar(i+1,j-1,k)-vstar(i+1,j,k)) ))/24.
     $       +u(i,j,k)*u(i,j,k)*tng(j)/radius*dzu(i,j,k)*areauu(j)
        ua(i,j,k)=ua(i,j,k)+uadv
        va(i,j,k)=va(i,j,k)+vadv

      enddo
      enddo
      enddo

#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 1,km
          write(*,*) 'adv',ua(i,j,k),va(i,j,k)
        enddo
        write(*,*)
      endif
#endif

      call wait_t2d_e1p(sfund,sfvnd,58)
      call wait_t2d_w1p(sfund,sfvnd,59)
      call wait_t3d_e1p(ub,vb,60,61)
      call wait_t3d_w1p(ub,vb,62,63)
#ifdef BISMFRIC
      call wait_t2d_e2p(sfund,sfvnd,64)
      call wait_t2d_w2p(sfund,sfvnd,65)
      call wait_t3d_e2p(ub,vb,66,67)
      call wait_t3d_w2p(ub,vb,68,69)
      call wait_t3d_e2(dzu,70)
      call wait_t3d_w2(dzu,71)

ccc   biharmonic Smagorinsky friction scheme
      do k = 1,km
      do j = 2,jml
      do i = 2,iml
         dux2(i,j,k)=ex(i,j,k)*ex(i-1,j,k)*dxr*csr(j)*
     &        (ub(i,j,k)-sfund(i,j)-ub(i-1,j,k)+sfund(i-1,j))
         dvx2(i,j,k)=ex(i,j,k)*ex(i-1,j,k)*dxr*csr(j)*
     &        (vb(i,j,k)-sfvnd(i,j)-vb(i-1,j,k)+sfvnd(i-1,j))
         duy2(i,j,k)=ex(i,j,k)*ex(i,j-1,k)*dyr*
     &        (ub(i,j,k)-sfund(i,j)-ub(i,j-1,k)+sfund(i,j-1))
         dvy2(i,j,k)=ex(i,j,k)*ex(i,j-1,k)*dyr*
     &        (vb(i,j,k)-sfvnd(i,j)-vb(i,j-1,k)+sfvnd(i,j-1))
         dux(i,j,k)=ex(i,j,k)*ex(i-1,j,k)*dxr*csr(j)*
     &        (ub(i,j,k)-ub(i-1,j,k))
         dvx(i,j,k)=ex(i,j,k)*ex(i-1,j,k)*dxr*csr(j)*
     &        (vb(i,j,k)-vb(i-1,j,k))
         duy(i,j,k)=ex(i,j,k)*ex(i,j-1,k)*dyr*
     &        (ub(i,j,k)-ub(i,j-1,k))
         dvy(i,j,k)=ex(i,j,k)*ex(i,j-1,k)*dyr*
     &        (vb(i,j,k)-vb(i,j-1,k))
      enddo
      enddo
      enddo

      do k = 1,km
      do j = 2,jml
      do i = 3,iml-1
         if(ex(i,j,k)==0 .and. ex(i+1,j,k)==1.)then
            dux(i,j,k)=dux(i+1,j,k)
            dvx(i,j,k)=dvx(i+1,j,k)
            dux2(i,j,k)=dux2(i+1,j,k)
            dvx2(i,j,k)=dvx2(i+1,j,k)
         elseif(ex(i,j,k)==0. .and. ex(i-1,j,k)==1.)then
            dux(i,j,k)=dux(i-1,j,k)
            dvx(i,j,k)=dvx(i-1,j,k)
            dux2(i,j,k)=dux2(i-1,j,k)
            dvx2(i,j,k)=dvx2(i-1,j,k)
         endif
      enddo
      enddo
      enddo

      do k = 1,km
      do j = 3,jml-1
      do i = 2,iml
         if(ex(i,j,k)==0. .and. ex(i,j+1,k)==1.)then
            duy(i,j,k)=duy(i,j+1,k)
            dvy(i,j,k)=dvy(i,j+1,k)
            duy2(i,j,k)=duy2(i,j+1,k)
            dvy2(i,j,k)=dvy2(i,j+1,k)
         elseif(ex(i,j,k)==0. .and. ex(i,j-1,k)==1.)then
            duy(i,j,k)=duy(i,j-1,k)
            dvy(i,j,k)=dvy(i,j-1,k)
            duy2(i,j,k)=duy2(i,j-1,k)
            dvy2(i,j,k)=dvy2(i,j-1,k)
         endif
      enddo
      enddo
      enddo

      call exch_t3d_n2p(duy,dvy,37,38)
      call exch_t3d_s1p(duy,dvy,39,40)
      call exch_t3d_e2p(dux,dvx,72,73)
      call exch_t3d_w1p(dux,dvx,74,75)
      call exch_t3d_n2p(duy2,dvy2,81,82)
      call exch_t3d_s1p(duy2,dvy2,83,84)
      call exch_t3d_e2p(dux2,dvx2,85,86)
      call exch_t3d_w1p(dux2,dvx2,87,88)

      do j = 2,jml-1
      do i = 2,iml-1
         bsmags_0(i,j)=0.
         bsmagn_0(i,j)=0.
         bsmagw_0(i,j)=0.
         bsmage_0(i,j)=0.
      enddo
      enddo

      call wait_t3d_n2p(duy,dvy,37,38)
      call wait_t3d_s1p(duy,dvy,39,40)
      call wait_t3d_e2p(dux,dvx,72,73)
      call wait_t3d_w1p(dux,dvx,74,75)
      call wait_t3d_n2p(duy2,dvy2,81,82)
      call wait_t3d_s1p(duy2,dvy2,83,84)
      call wait_t3d_e2p(dux2,dvx2,85,86)
      call wait_t3d_w1p(dux2,dvx2,87,88)

      do k = 1,km
      do j = 2,jml-1
      do i = 2,iml-1
         defrate_sw=dsqrt( (dux(i,j,k)-dvy(i,j,k))**2
     &        +(duy(i,j,k)+dvx(i,j,k))**2 )
         defrate_nw=dsqrt( (dux(i,j,k)-dvy(i,j+1,k))**2
     &        +(duy(i,j+1,k)+dvx(i,j,k))**2 )
         defrate_se=dsqrt( (dux(i+1,j,k)-dvy(i,j,k))**2
     &        +(duy(i,j,k)+dvx(i+1,j,k))**2 )
         defrate_ne=dsqrt( (dux(i+1,j,k)-dvy(i,j+1,k))**2
     &        +(duy(i,j+1,k)+dvx(i+1,j,k))**2 )

         bsmags(i,j,k)=(csmag)*delx(j)*delx(j)/pi/4.
     &    * dsqrt(defrate_sw+defrate_se)
         bsmagn(i,j,k)=(csmag)*delx(j)*delx(j)/pi/4.
     &    * dsqrt(defrate_nw+defrate_ne)
         bsmagw(i,j,k)=(csmag)*delx(j)*delx(j)/pi/2.
     &  * dsqrt((anhf(j)*defrate_sw+ashf(j)*defrate_nw)*areaur(j))
         bsmage(i,j,k)=(csmag)*delx(j)*delx(j)/pi/2.
     &  * dsqrt((anhf(j)*defrate_se+ashf(j)*defrate_ne)*areaur(j))

!        bsmags_0(i,j)=bsmags(i,j,k)*dzu(i,j,k)/(hrr(i,j)+1.-tex(i,j,k))
!        bsmagn_0(i,j)=bsmagn(i,j,k)*dzu(i,j,k)/(hrr(i,j)+1.-tex(i,j,k))
!        bsmagw_0(i,j)=bsmagw(i,j,k)*dzu(i,j,k)/(hrr(i,j)+1.-tex(i,j,k))
!        bsmage_0(i,j)=bsmage(i,j,k)*dzu(i,j,k)/(hrr(i,j)+1.-tex(i,j,k))
        bsmags_0(i,j)=bsmags_0(i,j)+
     &        bsmags(i,j,k)**2*dzu(i,j,k)/(hrr(i,j)+1.-ex(i,j,1))
        bsmagn_0(i,j)=bsmagn_0(i,j)+
     &        bsmagn(i,j,k)**2*dzu(i,j,k)/(hrr(i,j)+1.-ex(i,j,1))
        bsmagw_0(i,j)=bsmagw_0(i,j)+
     &        bsmagw(i,j,k)**2*dzu(i,j,k)/(hrr(i,j)+1.-ex(i,j,1))
        bsmage_0(i,j)=bsmage_0(i,j)+
     &        bsmage(i,j,k)**2*dzu(i,j,k)/(hrr(i,j)+1.-ex(i,j,1))

         d2ud2x(i,j,k)=dxr*csr(j)*ex(i,j,k)
     &        *(bsmage(i,j,k)*dux2(i+1,j,k)
     &        -bsmagw(i,j,k)*dux2(i,j,k))
         d2ud2y(i,j,k)=dyr*ex(i,j,k)
     &        *(bsmagn(i,j,k)*duy2(i,j+1,k)
     &        -bsmags(i,j,k)*duy2(i,j,k))
         d2vd2x(i,j,k)=dxr*csr(j)*ex(i,j,k)
     &        *(bsmage(i,j,k)*dvx2(i+1,j,k)
     &        -bsmagw(i,j,k)*dvx2(i,j,k))
         d2vd2y(i,j,k)=dyr*ex(i,j,k)
     &        *(bsmagn(i,j,k)*dvy2(i,j+1,k)
     &        -bsmags(i,j,k)*dvy2(i,j,k))

      enddo
      enddo
      enddo
c
      do j = 2,jml-1
      do i = 2,iml-1
         bsmags_0(i,j)=dsqrt(bsmags_0(i,j))
         bsmagn_0(i,j)=dsqrt(bsmagn_0(i,j))
         bsmagw_0(i,j)=dsqrt(bsmagw_0(i,j))
         bsmage_0(i,j)=dsqrt(bsmage_0(i,j))
      enddo
      enddo
c
      do k = 1,km
      do j = 3,jml-1
      do i = 3,iml-1
         dux(i,j,k)=dzumin(i,j,k)*dxr*csr(j)
     &        *(d2ud2x(i,j,k)-d2ud2x(i-1,j,k))
         dvx(i,j,k)=dzumin(i,j,k)*dxr*csr(j)
     &        *(d2vd2x(i,j,k)-d2vd2x(i-1,j,k))
         duy(i,j,k)=dzvmin(i,j,k)*dyr
     &        *(d2ud2y(i,j,k)-d2ud2y(i,j-1,k))
         dvy(i,j,k)=dzvmin(i,j,k)*dyr
     &        *(d2vd2y(i,j,k)-d2vd2y(i,j-1,k))
      enddo
      enddo
      enddo
#endif
c
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
#ifdef BISMFRIC
ccc   biharmonic friction
         uvis_b=-dy*( bsmage(i,j,k)*dux(i+1,j,k)
     &        -bsmagw(i,j,k)*dux(i,j,k) )
     &        -dx*cs(j)*( bsmagn(i,j,k)*duy(i,j+1,k)
     &        -bsmags(i,j,k)*duy(i,j,k) )
         vvis_b=-dy*( bsmage(i,j,k)*dvx(i+1,j,k)
     &        -bsmagw(i,j,k)*dvx(i,j,k) )
     &        -dx*cs(j)*( bsmagn(i,j,k)*dvy(i,j+1,k)
     &        -bsmags(i,j,k)*dvy(i,j,k) )
#else
         uvis_b=0.
         vvis_b=0.
#endif
c
ccc   laplacian friction
         uvis_l=fclfr*hduv*(
     $        dydxr*csr(j)*(
     &        dzumin(i,j,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i-1,j)+mask_lvis(i,j))
#endif
     &        *(ub(i-1,j,k)-ub(i,j,k)-sfund(i-1,j)+sfund(i,j))
     &        -dzumin(i+1,j,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i,j)+mask_lvis(i+1,j))
#endif
     &        *(ub(i,j,k)-ub(i+1,j,k)-sfund(i,j)+sfund(i+1,j))
     $        -( (dzu(i,j,k)-dzu(i-1,j,k))*xind(i,j,k)
     $        +(dzu(i,j,k)-dzu(i+1,j,k))*(1.-xind(i+1,j,k)) )
#ifdef MASKLVIS
     $        *mask_lvis(i,j)
#endif
     $        *2.d0*(ub(i,j,k)-sfund(i,j)) )
     $        +dxdyr*(
     &        dzvmin(i,j,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i,j-1)+mask_lvis(i,j))
#endif
     &        *(ub(i,j-1,k)-ub(i,j,k)-sfund(i,j-1)+sfund(i,j))
     &        -dzvmin(i,j+1,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i,j)+mask_lvis(i,j+1))
#endif
     &        *(ub(i,j,k)-ub(i,j+1,k)-sfund(i,j)+sfund(i,j+1))
     $        -( (dzu(i,j,k)-dzu(i,j-1,k))*yind(i,j,k)*cst(j)
     $        +(dzu(i,j,k)-dzu(i,j+1,k))*(1.-yind(i,j+1,k))*cst(j+1) )
#ifdef MASKLVIS
     $        *mask_lvis(i,j)
#endif
     $        *2.d0*(ub(i,j,k)-sfund(i,j)) ) )
        vvis_l=fclfr*hduv*(
     $        dydxr*csr(j)*(
     &        dzumin(i,j,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i-1,j)+mask_lvis(i,j))
#endif
     &        *(vb(i-1,j,k)-vb(i,j,k)-sfvnd(i-1,j)+sfvnd(i,j))
     &        -dzumin(i+1,j,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i,j)+mask_lvis(i+1,j))
#endif
     &        *(vb(i,j,k)-vb(i+1,j,k)-sfvnd(i,j)+sfvnd(i+1,j))
     $        -( (dzu(i,j,k)-dzu(i-1,j,k))*xind(i,j,k)
     $        +(dzu(i,j,k)-dzu(i+1,j,k))*(1.-xind(i+1,j,k)) )
     $        *2.d0*(vb(i,j,k)-sfvnd(i,j)) )
     $        +dxdyr*(
     &        dzvmin(i,j,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i,j-1)+mask_lvis(i,j))
#endif
     &        *(vb(i,j-1,k)-vb(i,j,k)-sfvnd(i,j-1)+sfvnd(i,j))
     &        -dzvmin(i,j+1,k)
#ifdef MASKLVIS
     &        *5.d-1*(mask_lvis(i,j)+mask_lvis(i,j+1))
#endif
     &        *(vb(i,j,k)-vb(i,j+1,k)-sfvnd(i,j)+sfvnd(i,j+1))
     $        -( (dzu(i,j,k)-dzu(i,j-1,k))*yind(i,j,k)*cst(j)
     $        +(dzu(i,j,k)-dzu(i,j+1,k))*(1.-yind(i,j+1,k))*cst(j+1) )
#ifdef MASKLVIS
     $        *mask_lvis(i,j)
#endif
     $        *2.d0*(vb(i,j,k)-sfvnd(i,j)) ) )
c
         ua(i,j,k)=ua(i,j,k)
     &        -dpdx(i,j,k)*dzu(i,j,k)*areauu(j)+uvis_b+uvis_l
         va(i,j,k)=va(i,j,k)
     &        -dpdy(i,j,k)*dzu(i,j,k)*areauu(j)+vvis_b+vvis_l

#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59
     &   .and. i==33 .and. j==51) then
          write(*,*) 'hor',ua(i,j,k),va(i,j,k),dpdx(i,j,k),dpdy(i,j,k),
     &     uvis_b,uvis_l,vvis_b,vvis_l
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
          write(*,*) 'hor',ua(i,j,k),va(i,j,k),dpdx(i,j,k),dpdy(i,j,k)
        enddo
        write(*,*)
      endif
#endif
      call wait_t3d_s1p(tvst,svst,4,5)
      call wait_t3d_w1p(tust,sust,76,77)
c
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
         ta(i,j,k)=ta(i,j,k)
     $        -tust(i,j,k)+tust(i-1,j,k)-tvst(i,j,k)+tvst(i,j-1,k)
         sa(i,j,k)=sa(i,j,k)
     $        -sust(i,j,k)+sust(i-1,j,k)-svst(i,j,k)+svst(i,j-1,k)
      enddo
      enddo
      enddo
#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 15,20
          write(*,*) 'hor-adv:',nkai,ta(i,j,k),
     &      tust(i,j,k),tust(i-1,j,k),tvst(i,j,k),tvst(i,j-1,k)
          write(*,*) '  ',sa(i,j,k),
     &      sust(i,j,k),sust(i-1,j,k),svst(i,j,k),svst(i,j-1,k)
        enddo
        write(*,*)
      endif
#endif
ccc   calculation of vertical mass flux (wl),
ccc   and vertical advection and diffusion of t and s
c
      do k = 1,km-1
      do j = 3,jml-1
      do i = 3,iml-1
         wl(i,j,k+1)=tex(i,j,k+1)*(wl(i,j,k)+0.5d0*(
     &         ustar(i,j,k)*(ex(i,j,k)+ex(i,j-1,k))
     &        -ustar(i-1,j,k)*(ex(i-1,j,k)+ex(i-1,j-1,k))
     &        +vstar(i,j,k)*(ex(i,j,k)+ex(i-1,j,k))
     &        -vstar(i,j-1,k)*(ex(i,j-1,k)+ex(i-1,j-1,k)) ))
      enddo
      enddo
      enddo
c
ccc   correction due to varible grid
c
      do k = 2,kadp
      do j = 3,jml-1
      do i = 3,iml-1
         wl(i,j,k) = wl(i,j,k)
     &        -( wl(i,j,1)+wflux(i,j)*areat(i,j,1) )
     &        *(1.d0-dep(k)/depadp)
      enddo
      enddo
      enddo

#ifdef DEBUG1
      if(nkai>482090.and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 1,km
          write(*,*) 'wl:',nkai,k,wl(i,j,k),
     &      ustar(i,j,k),ustar(i-1,j,k),vstar(i,j,k),vstar(i,j-1,k)
        enddo
        write(*,*)
      endif
#endif


ccc   mixed-layer scheme
      call noh_mix
c
ccc   isopycnal diffusion with gm scheme
      call gm_dif
#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 15,20
          write(*,*) 'gm:',nkai,ta(i,j,k),sa(i,j,k)
        enddo
        write(*,*)
      endif
#endif

      do j = 3,jml-2
      do k = 1,km-1
      do i = 3,iml-2
      if(tex(i,j,k+1).eq.1.)then
c
ccc   1st order
c        if(wl(i,j,k+1).gt.0.d0) then
c          facup = (1.d0-vupp)*0.5d0
c          facdn = (1.d0+vupp)*0.5d0
c        else
c          facup = (1.d0+vupp)*0.5d0
c          facdn = (1.d0-vupp)*0.5d0
c        endif
cc
c          tvt1=wl(i,j,k+1)*(facup*t(i,j,k)+facdn*t(i,j,k+1))
c          svt1=wl(i,j,k+1)*(facup*s(i,j,k)+facdn*s(i,j,k+1))
c
ccc   3rd order scheme
c
         tuu = tex(i,j,k-1)*tb(i,j,k-1)
     &        +(1.-tex(i,j,k-1))*tb(i,j,k)
         tdd = tex(i,j,k+2)*tb(i,j,k+2)
     &        +(1.-tex(i,j,k+2))*tb(i,j,k+1)
         suu = tex(i,j,k-1)*sb(i,j,k-1)
     &        +(1.-tex(i,j,k-1))*sb(i,j,k)
         sdd = tex(i,j,k+2)*sb(i,j,k+2)
     &        +(1.-tex(i,j,k+2))*sb(i,j,k+1)
c
        if(wl(i,j,k+1).gt.0.d0) then
          tc2 = 4.*(
     &     (dzt(i,j,k+1)+dzt(i,j,k+2))*(tb(i,j,k)-tb(i,j,k+1))
     &    -(dzt(i,j,k)+dzt(i,j,k+1))*(tb(i,j,k+1)-tdd)
c     &    +(dzt(i,j,k)+dzt(i,j,k+1))*(tb(i,j,k+1)-tdd)
     &    ) / ( 1.-tex(i,j,k+1)
     &   +(dzt(i,j,k+1)+dzt(i,j,k+2))*(dzt(i,j,k)+dzt(i,j,k+1))
     &   *(dzt(i,j,k)+2.*dzt(i,j,k+1)+dzt(i,j,k+2)) )
c
          sc2 = 4.*(
     &     (dzt(i,j,k+1)+dzt(i,j,k+2))*(sb(i,j,k)-sb(i,j,k+1))
     &    -(dzt(i,j,k)+dzt(i,j,k+1))*(sb(i,j,k+1)-sdd)
     &    ) / ( 1.-tex(i,j,k+1)
     &   +(dzt(i,j,k+1)+dzt(i,j,k+2))*(dzt(i,j,k)+dzt(i,j,k+1))
     &   *(dzt(i,j,k)+2.*dzt(i,j,k+1)+dzt(i,j,k+2)) )
        else
          tc2 = 4.*(
     &     (dzt(i,j,k)+dzt(i,j,k+1))*(tuu-tb(i,j,k))
     &    -(dzt(i,j,k-1)+dzt(i,j,k))*(tb(i,j,k)-tb(i,j,k+1))
     &    ) / (1.-tex(i,j,k)
     &   +(dzt(i,j,k)+dzt(i,j,k+1))*(dzt(i,j,k-1)+dzt(i,j,k))
     &   *(dzt(i,j,k-1)+2.*dzt(i,j,k)+dzt(i,j,k+1)))
          sc2 = 4.*(
     &     (dzt(i,j,k)+dzt(i,j,k+1))*(suu-sb(i,j,k))
     &    -(dzt(i,j,k-1)+dzt(i,j,k))*(sb(i,j,k)-sb(i,j,k+1))
     &    ) / (1.-tex(i,j,k)
     &   +(dzt(i,j,k)+dzt(i,j,k+1))*(dzt(i,j,k-1)+dzt(i,j,k))
     &   *(dzt(i,j,k-1)+2.*dzt(i,j,k)+dzt(i,j,k+1)))
        endif
c
c       tc2 = 0.
c       sc2 = 0.
c
        tc0 = (dzt(i,j,k+1)*tb(i,j,k)+dzt(i,j,k)*tb(i,j,k+1))
     &     /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k+1))
     &    -0.25*dzt(i,j,k)*dzt(i,j,k+1)*tc2
        sc0 = (dzt(i,j,k+1)*sb(i,j,k)+dzt(i,j,k)*sb(i,j,k+1))
     &     /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k+1))
     &    -0.25*dzt(i,j,k)*dzt(i,j,k+1)*sc2
        tc1 = 2.*(tb(i,j,k)-tb(i,j,k+1))
     &      /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k+1))
     &     -0.5*(dzt(i,j,k)-dzt(i,j,k+1))*tc2
        sc1 = 2.*(sb(i,j,k)-sb(i,j,k+1))
     &      /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k+1))
     &     -0.5*(dzt(i,j,k)-dzt(i,j,k+1))*sc2
c
ccc   in case of quick
c
c        tvt1 = wl(i,j,k+1)*tc0
c        svt1 = wl(i,j,k+1)*sc0
c
ccc   in case of quickest
c
       wc0 = wl(i,j,k+1)*c2dtts/(areat(i,j,k+1)+1.-tex(i,j,k+1))
c
       tvt1 = wl(i,j,k+1)*(tc0-0.5*wc0*tc1
     &    +(wc0*wc0-0.25*dzt(i,j,k)*dzt(i,j,k+1))/3.*tc2 )
       svt1 = wl(i,j,k+1)*(sc0-0.5*wc0*sc1
     &    +(wc0*wc0-0.25*dzt(i,j,k)*dzt(i,j,k+1))/3.*sc2 )
c
c
ccc   normal diffusion > see gm-dif
c
c        vdf=vdts(i,j,k+1)*areat(i,j,k+1)*2.d0
c     $    /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
c     &      *tex(i,j,k+1)
c          tvt2=-(tb(i,j,k)-tb(i,j,k+1))*vdf
c          svt2=-(sb(i,j,k)-sb(i,j,k+1))*vdf
cc
c        ta(i,j,k)=ta(i,j,k)+tvt2
c        ta(i,j,k+1)=ta(i,j,k+1)-tvt2
c        sa(i,j,k)=sa(i,j,k)+svt2
c        sa(i,j,k+1)=sa(i,j,k+1)-svt2
c
ccc   vertical advection at bottom is upstream
       if(tex(i,j,k+2).eq.0.)then
          if(wl(i,j,k+1).gt.0.d0) then
             facup = 0.d0
             facdn = 1.d0
          else
             facup = 1.d0
             facdn = 0.d0
          endif
          tvt1=wl(i,j,k+1)*(facup*t(i,j,k)+facdn*t(i,j,k+1))
          svt1=wl(i,j,k+1)*(facup*s(i,j,k)+facdn*s(i,j,k+1))
       endif
c
        ta(i,j,k)=ta(i,j,k)+tvt1
        ta(i,j,k+1)=ta(i,j,k+1)-tvt1
        sa(i,j,k)=sa(i,j,k)+svt1
        sa(i,j,k+1)=sa(i,j,k+1)-svt1
      endif
      enddo
      enddo
      enddo

!   short wave penetration
      do j = 3,jml-2
      do k = 1,kshrd
      do i = 3,iml-2
         if(k.le.kadp) then
            dpth = dep(k+1)*(1.+hcl(i,j)/depadp)
         else
            dpth = dep(k+1) + hcl(i,j)
         endif
#ifdef DEBUG
      if(dpth/35.<-1000.) then
        write(*,*) ip,i,j,k,dpth,hcl(i,j),dep(k+1)
      endif
#endif
         swpn = areat(i,j,k+1)*0.24*ddmnar*swr_o(i,j)
     &        *(0.58d0*exp(-dpth/35.d0)+0.42d0*exp(-dpth/2300.d0))

         ta(i,j,k)=ta(i,j,k)-swpn
         ta(i,j,k+1)=ta(i,j,k+1)+swpn

      enddo
      enddo
      enddo
#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 15,20
          write(*,*) 'rad:',nkai,ta(i,j,k)
        enddo
        write(*,*)
      endif
#endif
ccc   volume change
      do k = 1,kadp
      do j = 3,jml-2
      do i = 3,iml-2
         fact1 = dz(k)/depadp
         wuv =ex(i,j,k)/areauu(j)*
     &        (ashf(j)*(wl(i,j,1)+wl(i+1,j,1)
     &        +wflux(i,j)*areat(i,j,1)+wflux(i+1,j)*areat(i+1,j,1))
     &        +anhf(j)*(wl(i,j+1,1)+wl(i+1,j+1,1)
     &        +wflux(i,j+1)*areat(i,j+1,1)
     &        +wflux(i+1,j+1)*areat(i+1,j+1,1)))
         ua(i,j,k) = ua(i,j,k)-ub(i,j,k)*wuv*fact1
         va(i,j,k) = va(i,j,k)-vb(i,j,k)*wuv*fact1
         ta(i,j,k) = ta(i,j,k)
     &        -tb(i,j,k)*(wl(i,j,1)+wflux(i,j)*areat(i,j,1))*fact1
         sa(i,j,k) = sa(i,j,k)
     &        -sb(i,j,k)*(wl(i,j,1)+wflux(i,j)*areat(i,j,1))*fact1
      enddo
      enddo
      enddo

ccc   vertical advection and diffusion of u and v
      do j = 3,jml-2
      do k = 1,km-1
      do i = 3,iml-2
      if(k.le.exn(i,j)-1) then
c
c        vduv(i,j,k+1) = 1.d0
c
c          uvt=ex(i,j,k+1)*(0.5d0*(u(i,j,k)+u(i,j,k+1))
c     $         *(wl(i,j,k+1)*rar(i,j,k)+wl(i+1,j,k+1)*rar(i+1,j,k)
c     $         +wl(i,j+1,k+1)*rar(i,j+1,k)
c     $         +wl(i+1,j+1,k+1)*rar(i+1,j+1,k))
c     $         -vduv(i,j,k+1)*dzzr(k+1)*2.d0*areauu(j)
c     $         /(1.d0+dzr(k+1)*dzu(i,j,k+1))*(ub(i,j,k)-ub(i,j,k+1)))
c          vvt=ex(i,j,k+1)*(0.5d0*(v(i,j,k)+v(i,j,k+1))
c     $         *(wl(i,j,k+1)*rar(i,j,k)+wl(i+1,j,k+1)*rar(i+1,j,k)
c     $         +wl(i,j+1,k+1)*rar(i,j+1,k)
c     $         +wl(i+1,j+1,k+1)*rar(i+1,j+1,k))
c     $         -vduv(i,j,k+1)*dzzr(k+1)*2.d0*areauu(j)
c     $         /(1.d0+dzr(k+1)*dzu(i,j,k+1))*(vb(i,j,k)-vb(i,j,k+1)))
         uvt=ex(i,j,k+1)*(0.5d0*(u(i,j,k)+u(i,j,k+1))
     $        *(wl(i,j,k+1)*rar(i,j,k)+wl(i+1,j,k+1)*rar(i+1,j,k)
     $        +wl(i,j+1,k+1)*rar(i,j+1,k)
     $        +wl(i+1,j+1,k+1)*rar(i+1,j+1,k))
     $        -vduv(i,j,k+1)*2.d0*areauu(j)
     $        /(dzu(i,j,k)+dzu(i,j,k+1)+1.-ex(i,j,k))
     $        *(ub(i,j,k)-ub(i,j,k+1)))
         vvt=ex(i,j,k+1)*(0.5d0*(v(i,j,k)+v(i,j,k+1))
     $        *(wl(i,j,k+1)*rar(i,j,k)+wl(i+1,j,k+1)*rar(i+1,j,k)
     $        +wl(i,j+1,k+1)*rar(i,j+1,k)
     $        +wl(i+1,j+1,k+1)*rar(i+1,j+1,k))
     $        -vduv(i,j,k+1)*2.d0*areauu(j)
     $        /(dzu(i,j,k)+dzu(i,j,k+1)+1.-ex(i,j,k))
     $        *(vb(i,j,k)-vb(i,j,k+1)))
         ua(i,j,k)=ua(i,j,k)+uvt
         ua(i,j,k+1)=ua(i,j,k+1)-uvt
         va(i,j,k)=va(i,j,k)+vvt
         va(i,j,k+1)=va(i,j,k+1)-vvt

      else if(k.eq.exn(i,j)) then
ccc   advection of u and v along sloping bottom
c
         ua(i,j,k)=ua(i,j,k)+
     $        (0.5d0*wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        *(ex(i-1,j,k+1)*(u(i-1,j,k+1)+u(i,j,k))
     $        +ex(i-1,j-1,k+1)*(u(i-1,j-1,k+1)+u(i,j,k))
     $        +ex(i,j-1,k+1)*(u(i,j-1,k+1)+u(i,j,k)))
     $        +0.5d0*wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1)
     $        *(ex(i,j-1,k+1)*(u(i,j-1,k+1)+u(i,j,k))
     $        +ex(i+1,j-1,k+1)*(u(i+1,j-1,k+1)+u(i,j,k))
     $        +ex(i+1,j,k+1)*(u(i+1,j,k+1)+u(i,j,k)))
     $        +0.5d0*wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1)
     $        *(ex(i,j+1,k+1)*(u(i,j+1,k+1)+u(i,j,k))
     $        +ex(i-1,j+1,k+1)*(u(i-1,j+1,k+1)+u(i,j,k))
     $        +ex(i-1,j,k+1)*(u(i-1,j,k+1)+u(i,j,k)))
     $        +0.5d0*wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1)
     $        *(ex(i+1,j,k+1)*(u(i+1,j,k+1)+u(i,j,k))
     $        +ex(i+1,j+1,k+1)*(u(i+1,j+1,k+1)+u(i,j,k))
     $        +ex(i,j+1,k+1)*(u(i,j+1,k+1)+u(i,j,k))))
         va(i,j,k)=va(i,j,k)+
     $        (0.5d0*wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        *(ex(i-1,j,k+1)*(v(i-1,j,k+1)+v(i,j,k))
     $        +ex(i-1,j-1,k+1)*(v(i-1,j-1,k+1)+v(i,j,k))
     $        +ex(i,j-1,k+1)*(v(i,j-1,k+1)+v(i,j,k)))
     $        +0.5d0*wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1)
     $        *(ex(i,j-1,k+1)*(v(i,j-1,k+1)+v(i,j,k))
     $        +ex(i+1,j-1,k+1)*(v(i+1,j-1,k+1)+v(i,j,k))
     $        +ex(i+1,j,k+1)*(v(i+1,j,k+1)+v(i,j,k)))
     $        +0.5d0*wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1)
     $        *(ex(i,j+1,k+1)*(v(i,j+1,k+1)+v(i,j,k))
     $        +ex(i-1,j+1,k+1)*(v(i-1,j+1,k+1)+v(i,j,k))
     $        +ex(i-1,j,k+1)*(v(i-1,j,k+1)+v(i,j,k)))
     $        +0.5d0*wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1)
     $        *(ex(i+1,j,k+1)*(v(i+1,j,k+1)+v(i,j,k))
     $        +ex(i+1,j+1,k+1)*(v(i+1,j+1,k+1)+v(i,j,k))
     $        +ex(i,j+1,k+1)*(v(i,j+1,k+1)+v(i,j,k))))
      endif
c
      if(k.eq.exn(i-1,j-1)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        *(u(i,j,k+1)+u(i-1,j-1,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        *(v(i,j,k+1)+v(i-1,j-1,k))
      endif
      if(k.eq.exn(i-1,j)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*(wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        +wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1))
     $        *(u(i,j,k+1)+u(i-1,j,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*(wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        +wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1))
     $        *(v(i,j,k+1)+v(i-1,j,k))
      endif
      if(k.eq.exn(i-1,j+1)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1)
     $        *(u(i,j,k+1)+u(i-1,j+1,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1)
     $        *(v(i,j,k+1)+v(i-1,j+1,k))
      endif
      if(k.eq.exn(i,j-1)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*(wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        +wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1))
     $        *(u(i,j,k+1)+u(i,j-1,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*(wl(i,j,k+1)*rar(i,j,k)*rar(i,j,k+1)
     $        +wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1))
     $        *(v(i,j,k+1)+v(i,j-1,k))
      endif
      if(k.eq.exn(i,j+1)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*(wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1)
     $        +wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1))
     $        *(u(i,j,k+1)+u(i,j+1,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*(wl(i,j+1,k+1)*rar(i,j+1,k)*rar(i,j+1,k+1)
     $        +wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1))
     $        *(v(i,j,k+1)+v(i,j+1,k))
      endif
      if(k.eq.exn(i+1,j-1)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1)
     $        *(u(i,j,k+1)+u(i+1,j-1,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1)
     $        *(v(i,j,k+1)+v(i+1,j-1,k))
      endif
      if(k.eq.exn(i+1,j)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*(wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1)
     $        +wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1))
     $        *(u(i,j,k+1)+u(i+1,j,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*(wl(i+1,j,k+1)*rar(i+1,j,k)*rar(i+1,j,k+1)
     $        +wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1))
     $        *(v(i,j,k+1)+v(i+1,j,k))
      endif
      if(k.eq.exn(i+1,j+1)) then
         ua(i,j,k+1)=ua(i,j,k+1)
     $        -0.5d0*wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1)
     $         *(u(i,j,k+1)+u(i+1,j+1,k))
         va(i,j,k+1)=va(i,j,k+1)
     $        -0.5d0*wl(i+1,j+1,k+1)*rar(i+1,j+1,k)*rar(i+1,j+1,k+1)
     $        *(v(i,j,k+1)+v(i+1,j+1,k))
      endif
c
      enddo
      enddo
      enddo

c
#ifdef NESTED
! south
      if(ip_y .eq. 0) then
      do k=1,km
      do i=3,iml-2
      do j=3,jndg+2
            ta(i,j+2,k) = ta(i,j+2,k) +
     &           tex(i,j+2,k)*chft*(tsbt(i,j,k)-t(i,j+2,k))
     &           *areat(i,j+2,k)*dzt(i,j+2,k)
            sa(i,j+2,k) = sa(i,j+2,k) +
     &           tex(i,j+2,k)*chft*(ssbt(i,j,k)-s(i,j+2,k))
     &           *areat(i,j+2,k)*dzt(i,j+2,k)
      enddo
      enddo
      enddo
! avoid double counting
        jsn = jndg+5
      else
        jsn = 3
      endif
!north
      if(ip_y .eq. jpe-1) then
      do k=1,km
      do i=3,iml-2
      do j=3,jndg+2
            ta(i,jml-1-j,k) = ta(i,jml-1-j,k) +
     &        tex(i,jml-1-j,k)*chft*(tnbt(i,j,k)-t(i,jml-1-j,k))
     &       *areat(i,jml-1-j,k)*dzt(i,jml-1-j,k)
            sa(i,jml-1-j,k) = sa(i,jml-1-j,k) +
     &           tex(i,jml-1-j,k)*chft*(snbt(i,j,k)-s(i,jml-1-j,k))
     &           *areat(i,jml-1-j,k)*dzt(i,jml-1-j,k)
      enddo
      enddo
      enddo
! avoid double counting
        jen = jml-jndg-4
      else
        jen = jml-2
      endif
! west
      if(ip_x .eq. 0) then
      do k=1,km
      do j=jsn,jen
      do i=3,indg+2
            ta(i+2,j,k) = ta(i+2,j,k) +
     &           tex(i+2,j,k)*chft*(twbt(i,j,k)-t(i+2,j,k))
     &           *areat(i+2,j,k)*dzt(i+2,j,k)
            sa(i+2,j,k) = sa(i+2,j,k) +
     &           tex(i+2,j,k)*chft*(swbt(i,j,k)-s(i+2,j,k))
     &           *areat(i+2,j,k)*dzt(i+2,j,k)
      enddo
      enddo
      enddo
      endif
!  east
      if(ip_x .eq. ipe-1) then
      do k=1,km
      do j=jsn,jen
      do i=3,indg+2
            ta(iml-1-i,j,k) = ta(iml-1-i,j,k) +
     &           tex(iml-1-i,j,k)*chft*(tebt(i,j,k)-t(iml-1-i,j,k))
     &           *areat(iml-1-i,j,k)*dzt(iml-1-i,j,k)
            sa(iml-1-i,j,k) = sa(iml-1-i,j,k) +
     &           tex(iml-1-i,j,k)*chft*(sebt(i,j,k)-s(iml-1-i,j,k))
     &           *areat(iml-1-i,j,k)*dzt(iml-1-i,j,k)
      enddo
      enddo
      enddo

      endif
#endif

ccc   t/s/u/v forecast is until here
ccc   solving tri-diag matrix inversion for implicit vertical mixing
c
ccc   form matrix
c
ccc   k = 1
      do j = 3,jml-2
      do i = 3,iml-2
         tric1(i,j,1) = -vdts(i,j,2)*areat(i,j,2)*2.d0
     &        /(dzt(i,j,1)+dzt(i,j,2)+1.-tex(i,j,1))*voltr(i,j,1)
         trib1(i,j,1) = tex(i,j,1)/c2dtts-tric1(i,j,1)
c
         tric1t(i,j,1) = -vddt(i,j,2)*areat(i,j,2)*2.d0
     &        /(dzt(i,j,1)+dzt(i,j,2)+1.-tex(i,j,1))*voltr(i,j,1)
         trib1t(i,j,1) = tex(i,j,1)/c2dtts-tric1t(i,j,1)
c
         tric1s(i,j,1) = -vdds(i,j,2)*areat(i,j,2)*2.d0
     &        /(dzt(i,j,1)+dzt(i,j,2)+1.-tex(i,j,1))*voltr(i,j,1)
         trib1s(i,j,1) = tex(i,j,1)/c2dtts-tric1s(i,j,1)
c
         tric2(i,j,1) = -vduv(i,j,2)*areauu(j)*2.d0
     &        /(dzu(i,j,1)+dzu(i,j,2)+1.-ex(i,j,1))*volur(i,j,1)
         trib2(i,j,1) = ex(i,j,1)/c2dtuv-tric2(i,j,1)
c
      enddo
      enddo
c
      do k = 2,km-1
      do j = 3,jml-2
      do i = 3,iml-2
         tria1(i,j,k) = -tex(i,j,k)*vdts(i,j,k)*areat(i,j,k)*2.d0
     &        /(dzt(i,j,k-1)+dzt(i,j,k)+1.-tex(i,j,k-1))*voltr(i,j,k)
         tric1(i,j,k)=-tex(i,j,k+1)*vdts(i,j,k+1)*areat(i,j,k+1)*2.d0
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))*voltr(i,j,k)
         trib1(i,j,k) = tex(i,j,k)/c2dtts-tria1(i,j,k)-tric1(i,j,k)
c
         tria1t(i,j,k) = -tex(i,j,k)*vddt(i,j,k)*areat(i,j,k)*2.d0
     &        /(dzt(i,j,k-1)+dzt(i,j,k)+1.-tex(i,j,k-1))*voltr(i,j,k)
         tric1t(i,j,k)=-tex(i,j,k+1)*vddt(i,j,k+1)*areat(i,j,k+1)*2.d0
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))*voltr(i,j,k)
         trib1t(i,j,k) = tex(i,j,k)/c2dtts-tria1t(i,j,k)-tric1t(i,j,k)
c
         tria1s(i,j,k) = -tex(i,j,k)*vdds(i,j,k)*areat(i,j,k)*2.d0
     &        /(dzt(i,j,k-1)+dzt(i,j,k)+1.-tex(i,j,k-1))*voltr(i,j,k)
         tric1s(i,j,k)=-tex(i,j,k+1)*vdds(i,j,k+1)*areat(i,j,k+1)*2.d0
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))*voltr(i,j,k)
         trib1s(i,j,k) = tex(i,j,k)/c2dtts-tria1s(i,j,k)-tric1s(i,j,k)
c
         tria2(i,j,k) = -ex(i,j,k)*vduv(i,j,k)*areauu(j)*2.d0
     &        /(dzu(i,j,k-1)+dzu(i,j,k)+1.-ex(i,j,k-1))*volur(i,j,k)
         tric2(i,j,k) = -ex(i,j,k+1)*vduv(i,j,k+1)*areauu(j)*2.d0
     &        /(dzu(i,j,k)+dzu(i,j,k+1)+1.-ex(i,j,k))*volur(i,j,k)
         trib2(i,j,k) = ex(i,j,k)/c2dtuv-tria2(i,j,k)-tric2(i,j,k)
      enddo
      enddo
      enddo
c
ccc   k = km
      do j = 3,jml-2
      do i = 3,iml-2
         tria1(i,j,km) = -tex(i,j,km)*vdts(i,j,km)*areat(i,j,km)*2.d0
     &     /(dzt(i,j,km-1)+dzt(i,j,km)+1.-tex(i,j,km-1))*voltr(i,j,km)
         trib1(i,j,km) = tex(i,j,km)/c2dtts-tria1(i,j,km)
c
         tria1t(i,j,km) = -tex(i,j,km)*vddt(i,j,km)*areat(i,j,km)*2.d0
     &     /(dzt(i,j,km-1)+dzt(i,j,km)+1.-tex(i,j,km-1))*voltr(i,j,km)
         trib1t(i,j,km) = tex(i,j,km)/c2dtts-tria1t(i,j,km)
c
         tria1s(i,j,km) = -tex(i,j,km)*vdds(i,j,km)*areat(i,j,km)*2.d0
     &     /(dzt(i,j,km-1)+dzt(i,j,km)+1.-tex(i,j,km-1))*voltr(i,j,km)
         trib1s(i,j,km) = tex(i,j,km)/c2dtts-tria1s(i,j,km)
c
         tria2(i,j,km) = -ex(i,j,km)*vduv(i,j,km)*areauu(j)*2.d0
     &     /(dzu(i,j,km-1)+dzu(i,j,km)+1.-ex(i,j,km-1))*volur(i,j,km)
         trib2(i,j,km) = ex(i,j,km)/c2dtts-tria2(i,j,km)
      enddo
      enddo
c
ccc   solve equation
c
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,1).eq.1.d0) then
         ta(i,j,1) = ta(i,j,1)*voltr(i,j,1)
#ifdef ICE
         t_frz=s(i,j,1)*(a1_frz+a2_frz*dsqrt(s(i,j,1))
     &        +a3_frz*s(i,j,1))
         if(tb(i,j,1)+ta(i,j,1)*c2dtts.lt.t_frz)then
            ta1toice(i,j)=-(t_frz-(tb(i,j,1)+ta(i,j,1)*c2dtts))/c2dtts
            ta(i,j,1)=ta(i,j,1)-ta1toice(i,j)
         else
            ta1toice(i,j)=0.
         endif
#endif
c
         sa(i,j,1) = sa(i,j,1)*voltr(i,j,1)
c
         tribet1(i,j) = trib1(i,j,1)
         tribet1t(i,j) = trib1t(i,j,1)
         tribet1s(i,j) = trib1s(i,j,1)
c        td(i,j,1) = ta(i,j,1)/tribet1(i,j)
c        sd(i,j,1) = sa(i,j,1)/tribet1(i,j)
         td(i,j,1) = ta(i,j,1)/tribet1t(i,j)
         sd(i,j,1) = sa(i,j,1)/tribet1s(i,j)
c
      endif
c
      if(ex(i,j,1).eq.1.d0) then
         ua(i,j,1)= (ua(i,j,1)*areaur(j)+cor(j)*v(i,j,1)*dzu(i,j,1))
     &        /(dzua(i,j,1)+1.d0-ex(i,j,1))
         va(i,j,1)= (va(i,j,1)*areaur(j)-cor(j)*u(i,j,1)*dzu(i,j,1))
     &        /(dzua(i,j,1)+1.d0-ex(i,j,1))
c
         tribet2(i,j) = trib2(i,j,1)
         ud(i,j,1) = ua(i,j,1)/tribet2(i,j)
         vd(i,j,1) = va(i,j,1)/tribet2(i,j)
      endif
      enddo
      enddo

      do k = 2,km
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,k).eq.1.d0) then
         ta(i,j,k) = ta(i,j,k)*voltr(i,j,k)
         sa(i,j,k) = sa(i,j,k)*voltr(i,j,k)
c
c        trigam1(i,j,k) = tric1(i,j,k-1)/tribet1(i,j)
c        tribet1(i,j) = trib1(i,j,k) - tria1(i,j,k)*trigam1(i,j,k)
c        td(i,j,k) = (ta(i,j,k)-tria1(i,j,k)*td(i,j,k-1))/tribet1(i,j)
c        sd(i,j,k) = (sa(i,j,k)-tria1(i,j,k)*sd(i,j,k-1))/tribet1(i,j)
         trigam1t(i,j,k) = tric1t(i,j,k-1)/tribet1t(i,j)
         trigam1s(i,j,k) = tric1s(i,j,k-1)/tribet1s(i,j)
         tribet1t(i,j) = trib1t(i,j,k) - tria1t(i,j,k)*trigam1t(i,j,k)
         tribet1s(i,j) = trib1s(i,j,k) - tria1s(i,j,k)*trigam1s(i,j,k)
         td(i,j,k) = (ta(i,j,k)-tria1t(i,j,k)*td(i,j,k-1))/tribet1t(i,j)
         sd(i,j,k) = (sa(i,j,k)-tria1s(i,j,k)*sd(i,j,k-1))/tribet1s(i,j)
      endif
c
      if(ex(i,j,k).eq.1.d0) then ! add Coriolis component
         ua(i,j,k)= (ua(i,j,k)*areaur(j)+cor(j)*v(i,j,k)*dzu(i,j,k))
     &        /(dzua(i,j,k)+1.d0-ex(i,j,k))
         va(i,j,k)= (va(i,j,k)*areaur(j)-cor(j)*u(i,j,k)*dzu(i,j,k))
     &        /(dzua(i,j,k)+1.d0-ex(i,j,k))
c
         trigam2(i,j,k) = tric2(i,j,k-1)/tribet2(i,j)
         tribet2(i,j) = trib2(i,j,k) - tria2(i,j,k)*trigam2(i,j,k)
         ud(i,j,k) = (ua(i,j,k)-tria2(i,j,k)*ud(i,j,k-1))/tribet2(i,j)
         vd(i,j,k) = (va(i,j,k)-tria2(i,j,k)*vd(i,j,k-1))/tribet2(i,j)
      endif
      enddo
      enddo
      enddo
c
      do k = km-1,1,-1
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,k+1).eq.1.d0) then
         td(i,j,k) = td(i,j,k) - trigam1t(i,j,k+1)*td(i,j,k+1)
         sd(i,j,k) = sd(i,j,k) - trigam1s(i,j,k+1)*sd(i,j,k+1)
      endif
c
      if(ex(i,j,k+1).eq.1.d0) then
         ud(i,j,k) = ud(i,j,k) - trigam2(i,j,k+1)*ud(i,j,k+1)
         vd(i,j,k) = vd(i,j,k) - trigam2(i,j,k+1)*vd(i,j,k+1)
      endif
      enddo
      enddo
      enddo
#ifdef BNDTIDE
      td(:,:,:) = 0.
      sd(:,:,:) = 0.
#endif
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
         ta(i,j,k)=tex(i,j,k)*(tb(i,j,k)+td(i,j,k))
         sa(i,j,k)=tex(i,j,k)*(sb(i,j,k)+sd(i,j,k))
c
         ua(i,j,k)=ex(i,j,k)*(ub(i,j,k)+ud(i,j,k))
         va(i,j,k)=ex(i,j,k)*(vb(i,j,k)+vd(i,j,k))
c
         sfu(i,j)=sfu(i,j) + ua(i,j,k)*dzua(i,j,k)
     &        /(hrr(i,j)+hclua(i,j)+(1.-ex(i,j,k)))
         sfv(i,j)=sfv(i,j) + va(i,j,k)*dzua(i,j,k)
     &        /(hrr(i,j)+hclua(i,j)+(1.-ex(i,j,k)))
c
         zu(i,j) = zu(i,j) + dzu(i,j,k)*
     &        (ud(i,j,k)/c2dtts-cor(j)*v(i,j,k))
         zv(i,j) = zv(i,j) + dzu(i,j,k)*
     &        (vd(i,j,k)/c2dtts+cor(j)*u(i,j,k))
      enddo
      enddo
      enddo
c
#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59) then
        i=33
        j=51
        do k = 15,20
          write(*,*) 'last:',nkai,ta(i,j,k),sa(i,j,k),
     &      td(i,j,k),sd(i,j,k)
        enddo
        write(*,*)
      endif
#endif

#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do k=1,km
      do j=1,jml
      do i=1,indg+2
        twbt(i,j,k)=(twbc(i,j,k,mb1)
     $       +(twbc(i,j,k,mb2)-twbc(i,j,k,mb1))*cbf)
        swbt(i,j,k)=(swbc(i,j,k,mb1)
     $       +(swbc(i,j,k,mb2)-swbc(i,j,k,mb1))*cbf)
      enddo
      enddo
      enddo
      do k=1,km
      do j=1,jml
      do i=1,2
         ta(i+2,j,k)=twbt(i,j,k)*tex(i+2,j,k)
         sa(i+2,j,k)=swbt(i,j,k)*tex(i+2,j,k)
      enddo
      enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do k=1,km
      do j=1,jml
      do i=1,indg+2
        tebt(i,j,k)=(tebc(i,j,k,mb1)
     &       +(tebc(i,j,k,mb2)-tebc(i,j,k,mb1))*cbf)
        sebt(i,j,k) = (sebc(i,j,k,mb1)
     &       +(sebc(i,j,k,mb2)-sebc(i,j,k,mb1))*cbf)
      enddo
      enddo
      enddo
      do k=1,km
      do j=1,jml
      do i=1,2
         ta(iml-1-i,j,k)=tebt(i,j,k)*tex(iml-1-i,j,k)
         sa(iml-1-i,j,k)=sebt(i,j,k)*tex(iml-1-i,j,k)
      enddo
      enddo
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do k=1,km
      do i=1,iml
      do j=1,jndg+2
        tsbt(i,j,k)=(tsbc(i,j,k,mb1)
     $       +(tsbc(i,j,k,mb2)-tsbc(i,j,k,mb1))*cbf)
        ssbt(i,j,k)=(ssbc(i,j,k,mb1)
     $       +(ssbc(i,j,k,mb2)-ssbc(i,j,k,mb1))*cbf)
      enddo
      enddo
      enddo

      do k=1,km
      do i=1,iml
      do j=1,2
         ta(i,j+2,k)=tsbt(i,j,k)*tex(i,j+2,k)
         sa(i,j+2,k)=ssbt(i,j,k)*tex(i,j+2,k)
      enddo
      enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do k=1,km
      do i=1,iml
      do j=1,jndg+2
        tnbt(i,j,k)=(tnbc(i,j,k,mb1)
     $       +(tnbc(i,j,k,mb2)-tnbc(i,j,k,mb1))*cbf)
        snbt(i,j,k)=(snbc(i,j,k,mb1)
     $       +(snbc(i,j,k,mb2)-snbc(i,j,k,mb1))*cbf)
      enddo
      enddo
      enddo

      do k=1,km
      do i=1,iml
      do j=1,2
         ta(i,jml-1-j,k)=tnbt(i,j,k)*tex(i,jml-1-j,k)
         sa(i,jml-1-j,k)=snbt(i,j,k)*tex(i,jml-1-j,k)
      enddo
      enddo
      enddo
      endif
#endif


#ifdef ICE
      do j=3,jml-2
      do i=3,iml-2
         ta(i,j,1)=ta(i,j,1)+ta1toice(i,j)*c2dtts*tex(i,j,1)
      enddo
      enddo

      call ice_fwd
#endif

      if(matsno.eq.1.or.matsn2.eq.1) then
         do j = 1,jml
         do i = 1,iml
            hcl(i,j) = hcla(i,j)
            hclu(i,j) = hclua(i,j)
         enddo
         enddo
#ifdef NESTED
! west
      if(ip_x .eq. 0) then
         do j=1,jml
         do i=1,indg+2
            hclwbt(i,j)=hclawbt(i,j)
            wtwbt(i,j)=wtawbt(i,j)
         enddo
         enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
         do j=1,jml
!         do i=iml-indg-1,iml
         do i=1,indg+2
            hclebt(i,j)=hclaebt(i,j)
            wtebt(i,j)=wtaebt(i,j)
         enddo
         enddo
      endif
! south
      if(ip_y .eq. 0) then
         do i=1,iml
         do j=1,jndg+2
            hclsbt(i,j)=hclasbt(i,j)
            wtsbt(i,j)=wtasbt(i,j)
         enddo
         enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
         do i=1,iml
         do j=1,jndg+2
            hclnbt(i,j)=hclanbt(i,j)
            wtnbt(i,j)=wtanbt(i,j)
         enddo
         enddo
      endif
#endif
      else
         do j = 1,jml
         do i = 1,iml
            hclb(i,j) = hcl(i,j)
            hcl(i,j) = hcla(i,j)
            hclu(i,j) = hclua(i,j)
         enddo
         enddo

		 
#ifdef NESTED
! west
      if(ip_x .eq. 0) then
         do j=1,jml
         do i=1,indg+2
            hclbwbt(i,j)=hclwbt(i,j)
            hclwbt(i,j)=hclawbt(i,j)
            wtbwbt(i,j)=wtwbt(i,j)
            wtwbt(i,j)=wtawbt(i,j)
         enddo
         enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
         do j=1,jml
         do i=1,indg+2
            hclbebt(i,j)=hclebt(i,j)
            hclebt(i,j)=hclaebt(i,j)
            wtbebt(i,j)=wtebt(i,j)
            wtebt(i,j)=wtaebt(i,j)
         enddo
         enddo
      endif
! south
      if(ip_y .eq. 0) then
         do i=1,iml
         do j=1,jndg+2
            hclbsbt(i,j)=hclsbt(i,j)
            hclsbt(i,j)=hclasbt(i,j)
            wtbsbt(i,j)=wtsbt(i,j)
            wtsbt(i,j)=wtasbt(i,j)
         enddo
         enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
         do i=1,iml
         do j=1,jndg+2
            hclbnbt(i,j)=hclnbt(i,j)
            hclnbt(i,j)=hclanbt(i,j)
            wtbnbt(i,j)=wtnbt(i,j)
            wtnbt(i,j)=wtanbt(i,j)
         enddo
         enddo
      endif
#endif

      endif

	  
      do k = 1,km
      do j = 1,jml
      do i = 1,iml
         ua(i,j,k)=ex(i,j,k)*(ua(i,j,k)-sfu(i,j))
         va(i,j,k)=ex(i,j,k)*(va(i,j,k)-sfv(i,j))
      enddo
      enddo
      enddo


#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do k = 1, km
      do j = 3, jml-2
      do i = 1, 2
         ta(i+2,j,k)=twbt(i,j,k)*tex(i+2,j,k)
         sa(i+2,j,k)=swbt(i,j,k)*tex(i+2,j,k)
         ua(i+2,j,k) = ex(i+2,j,k)*( uwbc(i,j,k,mb1)
     &        + (uwbc(i,j,k,mb2) - uwbc(i,j,k,mb1)) * cbf )*factu
         va(i+2,j,k) = ex(i+2,j,k)*( vwbc(i,j,k,mb1)
     &        + (vwbc(i,j,k,mb2) - vwbc(i,j,k,mb1)) * cbf )*factu
#ifdef BNDTIDE
         ua(i+2,j,k) = ua(5,j,k)
         va(i+2,j,k) = va(5,j,k)
#endif
         zu(i+2,j) = 0.
         zv(i+2,j) = 0.
      enddo
      enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do k = 1, km
      do j = 1, jml
      do i = 1, 2
         ta(iml-1-i,j,k)=tebt(i,j,k)*tex(iml-1-i,j,k)
         sa(iml-1-i,j,k)=sebt(i,j,k)*tex(iml-1-i,j,k)
         ua(iml-2-i,j,k) = ex(iml-2-i,j,k)*( uebc(i,j,k,mb1)
     &        + (uebc(i,j,k,mb2) - uebc(i,j,k,mb1)) * cbf )*factu
         va(iml-2-i,j,k) = ex(iml-2-i,j,k)*( vebc(i,j,k,mb1)
     &        + (vebc(i,j,k,mb2) - vebc(i,j,k,mb1)) * cbf )*factu
#ifdef BNDTIDE
         ua(iml-2-i,j,k) = ua(iml-5,j,k)
         va(iml-2-i,j,k) = va(iml-5,j,k)
#endif
         zu(iml-2-i,j) = 0.
         zv(iml-2-i,j) = 0.
      enddo
        ua(iml-2,j,k) = 0.
        va(iml-2,j,k) = 0.
        zu(iml-2,j) = 0.
        zv(iml-2,j) = 0.
      enddo
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do k = 1, km
      do i = 1, iml
      do j = 1, 2
         ta(i,j+2,k)=tsbt(i,j,k)*tex(i,j+2,k)
         sa(i,j+2,k)=ssbt(i,j,k)*tex(i,j+2,k)
         ua(i,j+2,k) = ex(i,j+2,k)*( usbc(i,j,k,mb1)
     &        + (usbc(i,j,k,mb2) - usbc(i,j,k,mb1)) * cbf )*factu
         va(i,j+2,k) = ex(i,j+2,k)*( vsbc(i,j,k,mb1)
     &        + (vsbc(i,j,k,mb2) - vsbc(i,j,k,mb1)) * cbf )*factu
         zu(i,j+2) = 0.
         zv(i,j+2) = 0.
#ifdef BNDTIDE
         ua(i,j+2,k) = ua(i,5,k)
         va(i,j+2,k) = va(i,5,k)
#endif
      enddo
      enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do k = 1, km
      do i = 1, iml
      do j = 1, 2
         ta(i,jml-1-j,k)=tnbt(i,j,k)*tex(i,jml-1-j,k)
         sa(i,jml-1-j,k)=snbt(i,j,k)*tex(i,jml-1-j,k)
         ua(i,jml-2-j,k) = ex(i,jml-2-j,k)*( unbc(i,j,k,mb1)
     &        + (unbc(i,j,k,mb2) - unbc(i,j,k,mb1)) * cbf )*factu
         va(i,jml-2-j,k) = ex(i,jml-2-j,k)*( vnbc(i,j,k,mb1)
     &        + (vnbc(i,j,k,mb2) - vnbc(i,j,k,mb1)) * cbf )*factu
         zu(i,jml-2-j) = 0.
         zv(i,jml-2-j) = 0.
#ifdef BNDTIDE
         ua(i,jml-2-j,k) = ua(i,jml-5,k)
         va(i,jml-2-j,k) = va(i,jml-5,k)
#endif
      enddo
         ua(i,jml-2,k) = 0.
         va(i,jml-2,k) = 0.
         zu(i,jml-2) = 0.
         zv(i,jml-2) = 0.
      enddo
      enddo
      endif
#endif


!   transport things to prepare next step

      if(matsno.eq.1.or.matsn2.eq.1) then
         do k = 1,km
         do j = 1,jml
         do i = 1,iml
            t(i,j,k)=ta(i,j,k)
            s(i,j,k)=sa(i,j,k)
            u(i,j,k)=ua(i,j,k)
            v(i,j,k)=va(i,j,k)
         enddo
         enddo
         enddo
      else
         do k = 1,km
         do j = 1,jml
         do i = 1,iml
            tb(i,j,k)=t(i,j,k)
            sb(i,j,k)=s(i,j,k)
            ub(i,j,k)=u(i,j,k)
            vb(i,j,k)=v(i,j,k)
            t(i,j,k)=ta(i,j,k)
            s(i,j,k)=sa(i,j,k)
            u(i,j,k)=ua(i,j,k)
            v(i,j,k)=va(i,j,k)
         enddo
         enddo
         enddo
      endif

      return
      end