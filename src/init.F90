ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine init                                               c
c
c   $Id: init.F90 12 2008-12-11 11:12:34Z ishikawa $
c                                                                   c
c     This subroutine sets up the initial conditions for the        c
c     present job. according to the index 'nfirst' read in          c
c     'rdjobp', initial homogeneous state with rest or the final    c
c     state of the previous job is selectes.                        c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init

      use param
      use mod_mpi
      use mod_time,only : nfirst

      implicit none

#include "common.h"
c
      real*8 stbup !,wlg(im,jmg) !2009/12/06
      integer i,j,k,imin,ihour,iday,month,iyear
      real*8 aday,ayear
c
ccc   initilize all
c
      do k = 0,km+1
      do j = 1,jm
      do i = 1,im
        u(i,j,k)=0.d0
        v(i,j,k)=0.d0
        ustar(i,j,k)=0.d0
        vstar(i,j,k)=0.d0
        wl(i,j,k)=0.d0
        t(i,j,k)=0.d0
        s(i,j,k)=0.d0
c
        tke(i,j,k)=0.d0
c
        vdts(i,j,k)=0.d0
        vdtke(i,j,k)=0.d0
        vduv(i,j,k)=0.d0
        rho(i,j,k)=0.d0
        rhoo(i,j,k)=0.d0
        rhodf(i,j,k)=0.d0
        rhouf(i,j,k)=0.d0
        rhodh(i,j,k)=0.d0
        rhouh(i,j,k)=0.d0
c
        ub(i,j,k)=0.d0
        vb(i,j,k)=0.d0
c
        eust(i,j,k)=0.d0
        evst(i,j,k)=0.d0
        tudf(i,j,k)=0.d0
        tvdf(i,j,k)=0.d0
        sudf(i,j,k)=0.d0
        svdf(i,j,k)=0.d0
        ud(i,j,k)=0.d0
        vd(i,j,k)=0.d0
        td(i,j,k)=0.d0
        sd(i,j,k)=0.d0
        tked(i,j,k)=0.d0
c
#ifdef BISMFRIC
        dux(i,j,k)=0.d0
        duy(i,j,k)=0.d0
        dvx(i,j,k)=0.d0
        dvy(i,j,k)=0.d0
        dux2(i,j,k)=0.d0
        duy2(i,j,k)=0.d0
        dvx2(i,j,k)=0.d0
        dvy2(i,j,k)=0.d0
#endif

      enddo
      enddo
      enddo
c
      do j = 1,jm
      do i = 1,im
        hcl(i,j)=0.d0
        um(i,j)=0.d0
        vm(i,j)=0.d0
        sfund(i,j) = 0.d0
        sfvnd(i,j) = 0.d0
c
        umstar(i,j)=0.d0
        vmstar(i,j)=0.d0
c
        hclu(i,j)=0.d0
        hclua(i,j)=0.d0
c
        wsx(i,j)=0.d0
        wsy(i,j)=0.d0
c
        sfunb(i,j)=0.d0
        sfvnb(i,j)=0.d0
        uma(i,j)=0.d0
        vma(i,j)=0.d0
        hta(i,j)=0.d0
        hu(i,j)=0.d0
        hcla(i,j)=0.d0
        sfun(i,j)=0.d0
        sfvn(i,j)=0.d0

#ifdef ICE

        aice(i,j)=0.d0
        volice(i,j)=0.d0
        uice(i,j)=0.d0
        vice(i,j)=0.d0
#endif
      enddo
      enddo

#ifdef NCIO
      ncnt_rs=0
      ncnt_2d=0 !***
      ncnt_en=0 !***
      ncnt_rm=0
#endif
#if (defined(EXP2003) || defined(EXP2006)) && defined(ROKKA)
      nsec_scf(:)=0
      nsec_ws(:)=0
#endif
c
      if(nfirst) then
ccc   start from initial conditions (state at rest)
c
      ahour=0.d0
c
#ifdef NESTED
      call nestini
#else
      do k = 0,km+1
      do j = 1,jml
      do i = 1,iml
#ifdef GL11M
c        t(i,j,k)=tini(i,j,k)*tex(i,j,k)
c        s(i,j,k)=sini(i,j,k)*tex(i,j,k)
         t(i,j,k) = 0.5*(tref12(i,j,k,0)+tref12(i,j,k,1))*tex(i,j,k)
         s(i,j,k) = 0.5*(sref12(i,j,k,0)+sref12(i,j,k,1))*tex(i,j,k)
#endif
#ifdef EQ100M
        t(i,j,k)=tini(i,j,k)*tex(i,j,k)
        s(i,j,k)=sini(i,j,k)*tex(i,j,k)
#endif
#ifdef PC68M
!        t(i,j,k) = 0.5d0
!     &       *dble(tref12(i,j,k,12)+tref12(i,j,k,1))*tex(i,j,k)
!        s(i,j,k) = 0.5d0
!     &       *dble(sref12(i,j,k,12)+sref12(i,j,k,1))*tex(i,j,k)
! start from July 1
        t(i,j,k) = 0.5d0
     &       *dble(tref12(i,j,k,6)+tref12(i,j,k,7))*tex(i,j,k)
        s(i,j,k) = 0.5d0
     &       *dble(sref12(i,j,k,6)+sref12(i,j,k,7))*tex(i,j,k)
#endif
      enddo
      enddo
      enddo
#endif
c
      do k = 0,km+1
      do j = 1,jml
      do i = 1,iml
         tke(i,j,k)=tkemin*tex(i,j,k)
         vdts(i,j,k)=vdtsmn
         vdtke(i,j,k)=vdtsmn
         vduv(i,j,k)=vduvmn*ex(i,j,k)
      enddo
      enddo
      enddo
c
      do j = 1,jml
      do i = 1,iml
ccc   tke precondition
         tke(i,j,1) = 10.d0*tex(i,j,1)
         tke(i,j,2) = 5.0d0*tex(i,j,1)
         tke(i,j,3) = 1.d0 *tex(i,j,1)
         depml(i,j)=10.d2
      enddo
      enddo
c
      pm(1)=0.d0
      pd(0)=0.d0
      do k=1,km
         pd(k)=(9.81d-4)*dp(k)*1.036d0
         pm(k+1)=pm(k)+dz(k)*(9.81d-4)*1.036d0
      enddo
      pd(km+1)=pd(km)

      else
ccc   restart from the last data obtained by previous job
ccc   for new restart data
c
#ifdef PC68M
      if(ip.eq.imaster) then
         read(restart) last,ahour,nsec
      endif
      call bcast_int1(last,1)
      call bcast_dble1(ahour,1)
      call bcast_int1(nsec,1)

      if(ip.eq.imaster) then
         read(restart) pd,pm,ddmna,dmn
      endif
      call bcast_dble(pd,km+2)
      call bcast_dble(pm,km+2)
      call bcast_dble(ddmna,1)
      call bcast_dble(dmn,km+2)

      call read_3d_2d_r8(u,restart)
      call read_3d_2d_r8(v,restart)
      call read_3d_2d_r8(t,restart)
      call read_3d_2d_r8(s,restart)
      call read_3d_2d_r8(tke,restart)
c
      call read_2d_r8(hcl,restart)
      call read_2d_r8(um,restart)
      call read_2d_r8(vm,restart)
c
      call read_2d_r8(aice,restart)
      call read_2d_r8(volice,restart)
      call read_2d_r8(uice,restart)
      call read_2d_r8(vice,restart)
#else
#ifdef NESTED
      call nestini
#else
      if(ip.eq.imaster) then
        read(restart) last,ahour
      endif
      call bcast_int(last,1)
      call bcast_dble(ahour,1)
c
      if(ip.eq.imaster) then
        read(restart) pd,pm,ddmna,dmn
      endif
      call bcast_dble(pd,km+2)
      call bcast_dble(pm,km+2)
      call bcast_dble(ddmna,1)
      call bcast_dble(dmn,km+2)
c
      call read_3d_r8(u,restart)
      call read_3d_r8(v,restart)
      call read_3d_r8(t,restart)
      call read_3d_r8(s,restart)
      call read_3d_r8(tke,restart)
c
      call read_2d_r8(hcl,restart)
      call read_2d_r8(um,restart)
      call read_2d_r8(vm,restart)
c
      call read_2d_r8(aice,restart)
      call read_2d_r8(volice,restart)
      call read_2d_r8(uice,restart)
      call read_2d_r8(vice,restart)
#endif
#endif
c
      pd(0)=0.
      pd(km+1)=pd(km)

      if( nday == 360 ) then
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-iyear*8640.d0)*60.d0+.01d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
        write(*,*)last,' step : start from continue file'
        write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &     ' hour=',ihour,' min=',imin
        write(*,*)' '
      endif

      if(nint(dble(last)*dtts).ne.nint(ahour)*60*60)then
      if(ip.eq.imaster)then
         write(*,*)' '
         write(*,*)'time step interval (dtts) is changed!'
         write(*,*)'  from ',
     &        dnint(ahour*60.d0*60.d0/dble(last)),'sec'
      endif
      last=nint(ahour*60.d0*60.d0/dtts)
      if(ip.eq.imaster)then
         write(*,*)'  to   ',
     &        dnint(ahour*60.d0*60.d0/dble(last)),'sec'
         write(*,*)'starting step changes to ',last
         write(*,*)' '
      endif
      endif

      endif

      endif
c
ccc   end_of_read statments
c
      do k = 0,km+1
      do j = 1,jm
      do i = 1,im
        ub(i,j,k)=u(i,j,k)
        vb(i,j,k)=v(i,j,k)
        tb(i,j,k)=t(i,j,k)
        sb(i,j,k)=s(i,j,k)
        tkeb(i,j,k)=tke(i,j,k)
      enddo
      enddo
      enddo
c
ccc   unesco_density
c
      do k = 1,km
      do j = 1,jml
      do i = 1,iml
      if(tex(i,j,k).eq.1.d0)then
c
        rhoo(i,j,k)= 999.842594d0 + 6.793952d-2*t(i,j,k)
     &       - 9.095290d-3*t(i,j,k)*t(i,j,k)
     $       + 1.001685d-4*t(i,j,k)**3
     &       - 1.120083d-6*t(i,j,k)**4 + 6.536332d-9*t(i,j,k)**5
     &       + s(i,j,k)*(0.824493d0 - 4.0899d-3*t(i,j,k)
     &       + 7.6438d-5*t(i,j,k)*t(i,j,k)
     &       - 8.2467d-7*t(i,j,k)**3 + 5.3875d-9*t(i,j,k)**4)
     &       + dsqrt(s(i,j,k)**3)*(-5.72466d-3 + 1.0227d-4*t(i,j,k)
     &       - 1.6546d-6*t(i,j,k)*t(i,j,k)) + 4.8314d-4*s(i,j,k)**2
c
        rhoo(i,j,k)=rhoo(i,j,k)*1.d-3-1.d0
c
ccc   compute rho(s,theta,p)
c
        stbup = 1.965933d4 + 1.444304d2*t(i,j,k)
     $    - 1.706103d0*t(i,j,k)*t(i,j,k)
     &    + 9.648704d-3*t(i,j,k)**3  - 4.190253d-5*t(i,j,k)**4
     &    + s(i,j,k)*(52.84855d0 - 3.101089d-1*t(i,j,k)
     &    + 6.283263d-3*t(i,j,k)*t(i,j,k) -5.084188d-5*t(i,j,k)**3)
     &    + dsqrt(s(i,j,k)**3)*(3.886640d-1 + 9.085835d-3*t(i,j,k)
     &    - 4.619924d-4*t(i,j,k)*t(i,j,k))
     &    + pd(k)*(3.186519d0 + 2.212276d-2*t(i,j,k)
     &    - 2.984642d-4*t(i,j,k)*t(i,j,k) + 1.956415d-6*t(i,j,k)**3)
     &    + pd(k)*s(i,j,k)*(6.704388d-3  -1.847318d-4*t(i,j,k)
     &    + 2.059331d-7*t(i,j,k)*t(i,j,k))
     $    + 1.480266d-4*pd(k)*dsqrt(s(i,j,k)**3)
     &    + pd(k)*pd(k)*(2.102898d-4 - 1.202016d-5*t(i,j,k)
     &    + 1.394680d-7*t(i,j,k)*t(i,j,k))
     $    + pd(k)*pd(k)*s(i,j,k)*(-2.040237d-6
     &    + 6.128773d-8*t(i,j,k) + 6.207323d-10*t(i,j,k)*t(i,j,k))
c
        rho(i,j,k) = (rhoo(i,j,k)+1.d0)/(1.d0 - pd(k)/stbup) -1.d0
c
      endif
      enddo
      enddo
      enddo


      do j=2,jml
      do i=2,iml
          umstar(i,j)=.5d0*dy*((2.d0-ex(i,j-1,1))*um(i,j)
     $      +(2.d0-ex(i,j,1))*um(i,j-1))
          vmstar(i,j)=.5d0*dx*((2.d0-ex(i-1,j,1))*vm(i,j)
     $      +(2.d0-ex(i,j,1))*vm(i-1,j))*cs(j)
      enddo
      enddo
c
      call exch_2d_s1(vmstar,1)
      call exch_2d_w1(umstar,2)
      call wait_2d_s1(vmstar,1)
      call wait_2d_w1(umstar,2)

      do j=2,jml
      do i=2,iml
        wl(i,j,1)=.5d0*(umstar(i-1,j)*(ex(i-1,j,1)+ex(i-1,j-1,1))
     1       -umstar(i,j)*(ex(i,j,1)+ex(i,j-1,1))
     2       +vmstar(i,j-1)*(ex(i,j-1,1)+ex(i-1,j-1,1))
     3       -vmstar(i,j)*(ex(i,j,1)+ex(i-1,j,1)))*tex(i,j,1)
      enddo
      enddo
c
c      do j=1,jml
c      do i=1,iml
      do j=1,jm
      do i=1,im
        hclb(i,j)=hcl(i,j)
        umb(i,j)=um(i,j)
        vmb(i,j)=vm(i,j)
      enddo
      enddo
c
      do j=1,jml-1
      do i=1,iml-1
        hclu(i,j)=ex(i,j,1)*(ashf(j)*(hcl(i,j)+hcl(i+1,j))
     &    +anhf(j)*(hcl(i,j+1)+hcl(i+1,j+1)))*areaur(j)
        sfun(i,j)=ex(i,j,1)*um(i,j)
     &     /(hrr(i,j)+hclu(i,j)+(1.d0-ex(i,j,1)))
        sfvn(i,j)=ex(i,j,1)*vm(i,j)
     &     /(hrr(i,j)+hclu(i,j)+(1.d0-ex(i,j,1)))
        sfund(i,j) = sfun(i,j)
        sfvnd(i,j) = sfvn(i,j)
      enddo
      enddo
c
#ifdef ICE
      do j = 1,jml
      do i = 1,iml
        aiceb(i,j)=aice(i,j)
        voliceb(i,j)=volice(i,j)
        t_top(i,j) = t(i,j,1)
      enddo
      enddo
#endif
c

      return
      end
c

