! utilities for nesting
!  $Id: nest.F90 12 2008-12-11 11:12:34Z ishikawa $
c==================================================================

      subroutine nestini

#ifdef NESTED
      use param
      use mod_mpi
      use mod_time,only : nfirst
      implicit none
#include "common.h"

      real(4) ::  r4sf(im0,jm0,nday)!,r4o3(im0,jm0,0:km+1),r4o2(im0,jm0)
      real(8) ::  r8o3(im0,jm0,0:km+1),r8o2(im0,jm0)
      real(8) ::  twtr,tetr,tstr,tntr,twflx,wlev
      real(8) ::  twtr0,tetr0,tstr0,tntr0
#ifdef TRPCTL
      real(8),parameter :: trc_west=-3.d11,trc_south=-7.d11
      real(8),parameter :: trc_east=trc_west+trc_south
!      integer,parameter :: jwlim_tc=77,islim_tc=37
#endif
      real(8) :: cls,cln,clw,cle
      real(8) ::  cls0,cln0,clw0,cle0
      real(8) :: coefx,coefy
      real(8) :: dummy(4),adum,h_ice
      integer :: idum,jen,ien,isn,jsn,i,j,k,n,ii,jj,m,l
      integer :: month,iday,ihour,imin,isec
      integer :: ndcount/0/,nunum/00/
      real(8) ::  bttu,bttv,stbup
      integer :: ionum

#ifdef ROKKA
#ifdef EXP2003
      real(8) :: inihour=95.d0*24.d0   ! NWNPAC start from exp0303
#endif
#ifdef EXP2006
!      real(8) :: inihour=290.d0*24.d0   ! NWNPAC start from exp0610
      real(8) :: inihour=234.d0*24.d0   ! NWNPAC start from exp0608
      integer  :: num_inidt=5
#endif
#ifdef EXP2007
      real(8) :: inihour=202.d0*24.d0   ! NWNPAC start from exp0707 (July/22/2007)
      integer  :: num_inidt=5
#endif
#endif

!    NOTE:
!      ubd,vbd,uwbt,vwbt,uebt,vebt,usbt,vsbt,unbt,vnbt are all
!              "baroclinic" velocities
!
!  flag for off-line calculation
!   0: on-line mode, 1: offline read mode, 2: write mode(on-line)
! unit 61 flux data, unit 62 field data
      integer,parameter :: iflag_off =0

!  side bounday transport correction
!     0: no use   /=0: use
      integer,parameter :: iflag_tc=0


!  nday of open boundary model (climatology condition)
#ifdef NWNPAC
      integer :: nday_obc=365
!      integer :: nday_obc=360
#else
      integer :: nday_obc=365
#endif
#ifdef JP68M
      nday_obc=nday
#endif

      if(obc_clim.ne.1) nday_obc=nday

! initilization

      nrct=1

      do n = 1,nrbd
      do k = 0,km+1

      do j = 1,jm
      do i = 1,2
        uwbc(i,j,k,n) = 0.
        uebc(i,j,k,n) = 0.
        vwbc(i,j,k,n) = 0.
        vebc(i,j,k,n) = 0.
      enddo
      enddo

      do j = 1,2
      do i = 1,im
        usbc(i,j,k,n) = 0.
        unbc(i,j,k,n) = 0.
        vsbc(i,j,k,n) = 0.
        vnbc(i,j,k,n) = 0.
      enddo
      enddo

      do j = 1,jm
      do i = 1,indg+2
        twbc(i,j,k,n) = 0.
        tebc(i,j,k,n) = 0.
        swbc(i,j,k,n) = 0.
        sebc(i,j,k,n) = 0.
        tkewbc(i,j,k,n) = 0.
        tkeebc(i,j,k,n) = 0.
      enddo
      enddo

      do j = 1,jndg+2
      do i = 1,im
        tsbc(i,j,k,n) = 0.
        tnbc(i,j,k,n) = 0.
        ssbc(i,j,k,n) = 0.
        snbc(i,j,k,n) = 0.
        tkesbc(i,j,k,n) = 0.
        tkenbc(i,j,k,n) = 0.
      enddo
      enddo

      enddo
      enddo
!-------------- original model topography ---------------------------

      if(ip.eq.imaster) then
        read(omtopo)
        read(omtopo) exnom
      endif
      call bcast_int(exnom,im0*jm0)

! redundant !all pe calculate same
      do k=0,km+1
      do j=1,jm0
      do i=1,im0
         exom(i,j,k)=0.d0
         texom(i,j,k)=0.d0
         areatom(i,j,k)=0.d0
      enddo
      enddo
      enddo

      do j = 1,jm0
      do i = 1,im0
        texnom(i,j) = 0.
      enddo
      enddo

      do j=2,jm0
      do i=2,im0
         texnom(i,j)=max0(exnom(i,j),exnom(i-1,j),
     $        exnom(i,j-1),exnom(i-1,j-1))
      enddo
      enddo

#ifdef JP44
!+++++++++++++++++ cyclic condition for original model  +++++++++++++++
      do j=2,jm0
         texnom(1,j)=texnom(im0-3,j)
      enddo
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

      do k=1,km
      do j=1,jm0
      do i=1,im0
         if(k.le.exnom(i,j)) exom(i,j,k)=1.d0
         if(k.le.texnom(i,j)) texom(i,j,k)=1.d0
      enddo
      enddo
      enddo

      dxom=dxdeg0*radius/radian
      dyom=dydeg0*radius/radian
      ddyom=2.d0*radius*dsin(.5d0*dydeg0/radian)
      tanfi4om=dtan(.25d0*dydeg0/radian)
      dxddyom=dxom*ddyom

      do j=1,jm0
         sineom(j)=dsin((slatu0+dble(j-1)*dydeg0)/radian)
         csom(j)=dcos((slatu0+dble(j-1)*dydeg0)/radian)
         tngom(j)=sineom(j)/csom(j)
         ashfom(j)=.25d0*(1.d0+tngom(j)*tanfi4om)*csom(j)*dxddyom
         anhfom(j)=.25d0*(1.d0-tngom(j)*tanfi4om)*csom(j)*dxddyom
      enddo

      do j = 2,jm0
      do k = 1,km
      do i = 2,im0
         areatom(i,j,k)=texom(i,j,k)*
     &        ( (exom(i-1,j-1,k)+exom(i,j-1,k))*anhfom(j-1)
     &     +(exom(i-1,j,k)+exom(i,j,k))*ashfom(j) )
      enddo
      enddo
      enddo

#ifdef JP44X
! only for boundary pe
!   south
      if(ip_y .eq. 0) then
      do i=3,iml-2
         ii=(lw(ip_x)+i-5+ism*inc)/inc
         if(exnom(ii,jsm-1).ne.exn(i,3)) then
            write(*,*) exnom(ii,jsm-1),exn(i,3)
            write(*,*) 'Geometry error i=',i,' ,j=3'
            write(*,*) ii,jsm-1
            stop
         endif
         if(exnom(ii,jsm).ne.exn(i,4)) then
            write(*,*) exnom(ii,jsm),exn(i,4)
            write(*,*) 'Geometry error i=',i,' ,j=4'
            write(*,*) ii,jsm
            stop
         endif
      enddo
      endif
!   north
      if(ip_y .eq. jpe-1) then
      do i=3,iml-2
         ii=(lw(ip_x)+i-5+ism*inc)/inc
         if(exnom(ii,jsm+jg-1).ne.exn(i,jml-4)) then
            write(*,*) 'Geometry error i=',i,' ,j=jm-4'
            stop
         endif
         if(exnom(ii,jsm+jg).ne.exn(i,jml-3)) then
            write(*,*) 'Geometry error i=',i,' ,j=jm-3'
            stop
         endif
      enddo
      endif

!   west
      if(ip_x .eq. 0) then
      do j=3,jml-2
         jj=(ls(ip_y)+j-5+jsm*jnc)/jnc
         if(exnom(ism-1,jj).ne.exn(3,j)) then
            write(*,*) 'Geometry error i=3, j=',j
            write(*,*) ism-1,jj
            stop
         endif
         if(exnom(ism,jj).ne.exn(4,j)) then
            write(*,*) exom(ism,jj,1),ex(4,j,1)
            write(*,*) 'Geometry error i=4, j=',j
            stop
         endif
      enddo
      endif
!   east
      if(ip_x .eq. ipe-1) then
      do j=3,jml-2
         jj=(ls(ip_y)+j-5+jsm*jnc)/jnc
         if (exnom(ism+ig-1,jj).ne.exn(iml-4,j)) then
            write(*,*) 'Geometry error i=im-4,  j=',j
            stop
         endif
         if (exnom(ism+ig,jj).ne.exn(iml-3,j)) then
            write(*,*) 'Geometry error i=im-3,  j=',j
            stop
         endif
      enddo
      endif
#endif

#ifdef TRPCTL
! Transport control through Tsushima strait

      do n = 1,4
        dummy(n) = 0.
      enddo
!  west
      if(ip_x .eq. 0) then
        do j=3,jml-2
           jj=ls(ip_y)-1+j
           if(jj<=jlim_tc) then
              dummy(1)=dummy(1)+ex(4,j,1)*dy
           endif
        enddo
        if(ip_y .eq. 0) then
          dummy(1)=dummy(1)+ex(4,2,1)*dy
        endif
      endif
!  east
      if(ip_x .eq. ipe-1) then
        do j=3,jml-2
          dummy(2)=dummy(2)+ex(iml-4,j,1)*dy
        enddo
        if(ip_y .eq. 0) then
          dummy(2)=dummy(2)+ex(iml-4,2,1)*dy
        endif
      endif
!  south
      if(ip_y .eq. 0) then
        do i=3,iml-2
           ii=lw(ip_x)-1+i
           if(ii<=ilim_tc) then
              dummy(3)=dummy(3)+ex(i,4,1)*dx*cs(4)
           endif
        enddo
        if(ip_x .eq. 0) then
          dummy(3)=dummy(3)+ex(2,4,1)*dx*cs(4)
        endif
      endif
!   north
      if(ip_y .eq. jpe-1) then
        do i=3,iml-2
          dummy(4)=dummy(4)+ex(i,jml-4,1)*dx*cs(jml-4)
        enddo
        if(ip_x .eq. 0) then
          dummy(4)=dummy(4)+ex(2,jml-4,1)*dx*cs(jml-4)
        endif
      endif

      call sum_all(dummy,4)
      clw=dummy(1)
      cle=dummy(2)
      cls=dummy(3)
      cln=dummy(4)

      umc_w=trc_west/clw
      umc_e=trc_east/cle
      vmc_s=trc_south/cls

      if(ip.eq.imaster) then
         write(*,*) 'transport control'
         write(*,*) 'um west',umc_w
         write(*,*) 'um east',umc_e
         write(*,*) 'vm south',vmc_s
      endif

! redundant
      clw0=0.d0
      cle0=0.d0
      cls0=0.d0
      cln0=0.d0
      do j=jsm,jsm+jg-1
         clw0=clw0+exom(ism,j,1)*dyom
         cle0=cle0+exom(ism+ig-1,j,1)*dyom
      enddo
      do i=ism,ism+ig-1
         cls0=cls0+exom(i,jsm,1)*dxom*cs(2)
         cln0=cln0+exom(i,jsm+jg-1,1)*dxom*cs(jm-2)
      enddo

#endif

      do k=0,km+1
      do n=1,ndpt
         ipt(n,k)=0
         jpt(n,k)=0
         ipu(n,k)=0
         jpu(n,k)=0
      enddo
      enddo

      do k=1,km
         isumt(k)=0
         isumu(k)=0
         do j=1,jm0
         do i=1,im0
            if(texom(i,j,k).ne.0.d0) then
               isumt(k)=isumt(k)+1
               ipt(isumt(k),k)=i
               jpt(isumt(k),k)=j
            endif
            if(exom(i,j,k).ne.0.d0) then
               isumu(k)=isumu(k)+1
               ipu(isumu(k),k)=i
               jpu(isumu(k),k)=j
            endif
         enddo
         enddo
      enddo

#ifdef JP44
!==================== surface bounday condition ======================
!------------------------------- wsx_d --------------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(wsx_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplws(wsx_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(ex(i,j,1).ne.0.d0.and.wsx_d(i,j,m).eq.999.0) then
            write(*,*) 'ip=',ip,'i=',i,' j=',j,' m=',m,
     &      ex(i,j,1),wsx_d(i,j,m)
            write(*,*) 'unexpected value for wsx_d'
            stop
         endif
#endif
         wsx_d(i,j,m)=sngl(wsx_d(i,j,m)*ex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(wsx_d,61)
      endif
      if(ip .eq. imaster) write(*,*) 'wsx_d is finished'
!----------------- wsy_d ------------------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(wsy_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplws(wsy_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(ex(i,j,1).ne.0.d0.and.wsy_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for wsy_d'
            stop
         endif
#endif
         wsy_d(i,j,m)=sngl(wsy_d(i,j,m)*ex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(wsy_d,61)
      endif
      if(ip .eq. imaster) write(*,*) 'wsy_d is finished'

!------------------------- sc_wind_d --------------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(sc_wind_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplsf(sc_wind_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(tex(i,j,1).ne.0.d0.and.sc_wind_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for sc_wind_d'
            stop
         endif
#endif
         sc_wind_d(i,j,m)=sngl(sc_wind_d(i,j,m)*tex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(sc_wind_d,61)
      endif
      if(ip .eq. imaster) write(*,*) 'sc_wind_d is finished'

!------------------- stdev_w_d ----------------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(stdev_w_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplsf(stdev_w_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(tex(i,j,1).ne.0.d0.and.stdev_w_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for stdev_w_d'
            stop
         endif
#endif
         stdev_w_d(i,j,m)=sngl(stdev_w_d(i,j,m)*tex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(stdev_w_d,61)
      endif
      if(ip .eq. imaster) write(*,*) 'stdev_w_d is finished'


!------------------------- dpt_t2m_d ----------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(dpt_t2m_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplsf(dpt_t2m_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(tex(i,j,1).ne.0.d0.and.dpt_t2m_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for dpt_t2m_d'
            stop
         endif
#endif
         dpt_t2m_d(i,j,m)=sngl(dpt_t2m_d(i,j,m)*tex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(dpt_t2m_d,61)
      endif
      if(ip .eq. imaster) write(*,*) 'dep_t2m_d is finished'

!----------------------- temp_2m_d ---------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(temp_2m_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplsf(temp_2m_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(tex(i,j,1).ne.0.d0.and.temp_2m_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for temp_2m_d'
            stop
         endif
#endif
         temp_2m_d(i,j,m)=sngl(temp_2m_d(i,j,m)*tex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(temp_2m_d,61)
      endif

      if(ip .eq. imaster) write(*,*) 'temp_2m_d is finished'
!--------------------------- cloud_d ---------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(cloud_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplsf(cloud_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(tex(i,j,1).ne.0.d0.and.cloud_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for cloud_d'
            stop
         endif
#endif
         cloud_d(i,j,m)=sngl(cloud_d(i,j,m)*tex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(cloud_d,61)
      endif
      if(ip .eq. imaster) write(*,*) 'cloud_d is finished'

!----------------------- tot_sol_d ----------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(tot_sol_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplsf(tot_sol_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(tex(i,j,1).ne.0.d0.and.tot_sol_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for tot_sol_d'
            stop
         endif
#endif
         tot_sol_d(i,j,m)=sngl(tot_sol_d(i,j,m)*tex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(tot_sol_d,61)
      endif
      if(ip .eq. imaster) write(*,*) 'tot_sol_d is finished'

!-------------------------- wflux_d ------------------------
      if(ip.eq.imaster) then
        read(bnddt) r4sf
      endif

      if(iflag_off .eq. 1) then
        call read_fd_nest(wflux_d,61)
      else

      call bcast_real(r4sf,im0*jm0*nday)
      call intplsf(wflux_d,r4sf)

      do m=1,nday
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG1
         if(tex(i,j,1).ne.0.d0.and.wflux_d(i,j,m).eq.999.0) then
            write(*,*) 'unexpected value for wflux_d'
            stop
         endif
#endif
         wflux_d(i,j,m)=sngl(wflux_d(i,j,m)*tex(i,j,1))
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_fd_nest(wflux_d,61)
      endif

      if(ip .eq. imaster) write(*,*) 'wflux_d is finished'

#endif
!==============initial and size boundary condition ==============

      if(ip.eq.imaster) then
        write(*,*) 'intepolation field variables'
      endif
      do 300 m=1,nrbd
#ifndef ASSIM
      if(ip.eq.imaster) then
        write(*,*) 'time ;',m
      endif
#else
      if(ip.eq.imaster .and. mod(m,100) .eq. 0) then
        write(*,*) 'time ;',m
      endif
#endif
!-----------------------------------------------------------

!      if(ip.eq.imaster) then
!        read(omunt) idum,ahourbd
!        read(omunt) pdom,pm,ddmnar
!        dummy(1) = ahourbd
!        dummy(2) = ddmnar
!        write(*,*) 'nkaibd,ahourbd',idum,ahourbd
!      endif
!      call bcast_int(idum,1)
!      nkaibd(m) = idum
!      dummy(1) = ahourbd
!      dummy(2) = ddmnar
!      call bcast_dble(dummy,2)
!      ahourbd = dummy(1)
!      ddmnar = dummy(2)

!---------------- output new version -----------------------
!      if(ip.eq.imaster) then
!        read(omunt) adum
!        read(omunt) pdom,ddmnar
!        dummy(1) = adum
!        dummy(2) = ddmnar
!      endif
!      call bcast_int(idum,1)
!      dummy(1) = adum
!      dummy(2) = ddmnar
!      call bcast_dble(dummy,2)
!      ahourbd(m) = dummy(1)
!      ddmnar = dummy(2)

!-----------------------------------------------------------
      if(ip.eq.imaster) then
!        read(omunt) idum,adum
!        read(omunt) pdom,pm,ddmnar
         do
            read(omunt,iostat=ionum) adum
            select case(ionum)
            case default
               write(*,*) 'error number=',ionum
               write(*,*) 'field data reading err : ahour'
               stop
            case(0)
               exit
            case(-1)
               write(*,*) 'next obc file'
               omunt=omunt+1
            end select
         enddo
         read(omunt) pdom,ddmnar
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : pd,ddmnar'
           stop
        endif
#endif
        dummy(1) = adum
        dummy(2) = ddmnar
#ifndef ASSIM
        write(*,*) 'nkaibd,ahourbd',idum,adum
#endif
      endif
!      call bcast_int1(idum,1)
!      nkaibd(m) = idum
!      dummy(1) = adum
!      dummy(2) = ddmnar
      call bcast_dble(dummy,2)
      ahourbd(m) = dummy(1)*dble(nday)/dble(nday_obc)
#ifdef ROKKA
      ahourbd(m)=ahourbd(m)+inihour
#endif
#if defined(EXP2003)
         ahourbd(m)=ahourbd(m)-dble(360*24*2)
#endif
!      ddmnar = dummy(2)

      if(m.gt.1 .and. ahourbd(m).lt.ahourbd(m-1)) then
         if(obc_clim.eq.1) then

         ahourbd(m)=ahourbd(m)+24.d0*dble(nday)
#ifdef DEBUG
         if(ahourbd(m).lt.ahourbd(m-1) .or.
     $        ahourbd(m).ge.ahourbd(m-1)+24.d0*dble(nday)) then
            if(ip.eq.imaster) write(*,*) 'side boundary ahour err'
            stop
         endif
#endif
         else
            if(ip.eq.imaster) write(*,*) 'side boundary ahour err'
            stop
         endif

      endif

      call bcast_dble(pdom,km+2)
      call bcast_dble(pm,km+2)

!------------------------ t ----------------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o3
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : t'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_3d_r8(tbd,62)
      else

      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplt(tbd,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
         tbd(i,j,k)=tbd(i,j,k)*tex(i,j,k)
#ifdef DEBUG
         if(nint(tbd(i,j,k))==999) then
            write(*,*) 'unexpected value for tbd at'
     $           ,i+lw(ip_x)-1,j+ls(ip_y)-1,k
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_3d8_r8(tbd,62)
      endif

!------------------------- s -------------------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o3
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : s'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_3d_r8(sbd,62)
      else


      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplt(sbd,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
         sbd(i,j,k)=sbd(i,j,k)*tex(i,j,k)
#ifdef DEBUG
         if(nint(sbd(i,j,k))==999) then
            write(*,*) 'unexpected value for sbd,'
     $           ,i+lw(ip_x)-1,j+ls(ip_y)-1,k
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_3d8_r8(sbd,62)
      endif

!-------------------------- tke -----------------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o3
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : tke'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_3d_r8(tkebd,62)
      else

      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplt(tkebd,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(tex(i,j,k)>0.5 .and. nint(tkebd(i,j,k))==999) then
            write(*,*) 'unexpected value for tkebd at'
     $           ,i+lw(ip_x)-1,j+ls(ip_y)-1,k
#ifndef DEBUG20
            stop
#endif
         endif
#endif
         tkebd(i,j,k)=dmin1(dmax1(tkebd(i,j,k),tkemin),1.d6)*tex(i,j,k)
      enddo
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_3d8_r8(tkebd,62)
      endif


!---------------------------- hcl -----------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o2
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : hcl'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_2d_r8(hclbd,62)
      else

      call bcast_dble(r8o2,im0*jm0)
      call intplh(hclbd,r8o2)

      do j = 1, jml
      do i = 1, iml
         hclbd(i,j)=hclbd(i,j)*tex(i,j,1)
#ifdef DEBUG
         if(nint(hclbd(i,j))==999) then
            write(*,*) 'unexpected value for hclbd at'
     $           ,i+lw(ip_x)-1,j+ls(ip_y)-1
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo
      endif


      if(iflag_off .eq. 2) then
        call write_2d8_r8(hclbd,62)
      endif

!      write(96,*) hclbd
!     read(omunt) !htom

!-------------------------- um ----------------------------------
      if(ip.eq.imaster) then
        read(omunt) umom
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : um'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_2d_r8(umbd,62)
      else

      call bcast_dble(umom,im0*jm0)
      call intplum(umbd,umom)

      do j = 1, jml
      do i = 1, iml
         umbd(i,j)=umbd(i,j)*ex(i,j,1)
#ifdef DEBUG
         if(umbd(i,j)==999.d0) then
            write(*,*) 'unexpected value for umbd at',
     $           i+lw(ip_x)-1,j+ls(ip_y)-1
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo
      endif


      if(iflag_off .eq. 2) then
        call write_2d8_r8(umbd,62)
      endif

!--------------------------- vm ---------------------------------
      if(ip.eq.imaster) then
        read(omunt) vmom
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : vm'
!           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_2d_r8(vmbd,62)
      else

      call bcast_dble(vmom,im0*jm0)
      call intplvm(vmbd,vmom)

      do j = 1, jml
      do i = 1, iml
         vmbd(i,j)=vmbd(i,j)*ex(i,j,1)
#ifdef DEBUG
         if(vmbd(i,j)==999.d0) then
            write(*,*) 'unexpected value for vmbd at',
     $           i+lw(ip_x)-1,j+ls(ip_y)-1
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_2d8_r8(vmbd,62)
      endif

!------------------- vertical velocity on k=1 --------------------
      if(iflag_off .eq. 1) then
        call read_2d_r8(wtbd,62)
      else

      do j=3,jm0-1
      do i=3,im0-1
         wtom(i,j)=texom(i,j,1)*.25d0*(
     &        dyom*((exom(i-1,j,1)+exom(i-1,j-1,1))
     $        *((2.d0-exom(i-1,j-1,1))*umom(i-1,j)
     $        +(2.d0-exom(i-1,j,1))*umom(i-1,j-1))
     $        -(exom(i,j,1)+exom(i,j-1,1))
     &        *((2.d0-exom(i,j-1,1))*umom(i,j)
     $        +(2.d0-exom(i,j,1))*umom(i,j-1)))
     $        +dxom*((exom(i,j-1,1)+exom(i-1,j-1,1))
     $        *((2.d0-exom(i-1,j-1,1))*vmom(i,j-1)
     $        +(2.d0-exom(i,j-1,1))*vmom(i-1,j-1))*csom(j-1)
     $        -(exom(i,j,1)+exom(i-1,j,1))
     $        *((2.d0-exom(i-1,j,1))*vmom(i,j)
     $        +(2.d0-exom(i,j,1))*vmom(i-1,j))*csom(j)) )
     $        /(areatom(i,j,1)+1.d0-texom(i,j,1))
      enddo
      enddo

#ifdef JP44
!++++++++++++++++ cyclic original model ++++++++++++++++++
      do j=3,jm0-1
         wtom(1,j)    = wtom(im0-3,j)
         wtom(2,j)    = wtom(im0-2,j)
         wtom(im0-1,j) = wtom(3,j)
         wtom(im0,j)   = wtom(4,j)
      enddo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

      call intplh(wtbd,wtom)

      do j = 1, jml
      do i = 1, iml
         wtbd(i,j)=wtbd(i,j)*tex(i,j,1)
#ifdef DEBUG
         if(wtbd(i,j)==999.d0) then
            write(*,*) 'unexpected value for wtbd'
            stop
         endif
#endif
      enddo
      enddo
      endif

      if(iflag_off .eq. 2) then
        call write_2d8_r8(wtbd,62)
      endif


!------------------- total volume change -------------------------

      if(iflag_off .ne. 1) then
      do n = 1,4
        dummy(n) = 0.
      enddo

!  west
      if(ip_x .eq. 0 ) then
        if(ip_y .eq. jpe-1) then
          jen = jml-4
        else
          jen = jml-2
        endif
        if(ip_y .eq. 0) then
          jsn = 4
        else
          jsn = 3
        endif
        do j=jsn,jen
          dummy(1)=dummy(1)+umbd(4,j)*dy
        enddo
      endif
!  east
      if(ip_x .eq. ipe-1) then
        if(ip_y .eq. jpe-1) then
          jen = jml-4
        else
          jen = jml-2
        endif
        if(ip_y .eq. 0) then
          jsn = 4
        else
          jsn = 3
        endif
        do j=jsn,jen
          dummy(2)=dummy(2)+umbd(iml-3,j)*dy
        enddo
        if(ip_y .eq. 0 ) then
          dummy(2)=dummy(2)+umbd(iml-3,2)*dy
        endif
      endif
!  south
      if(ip_y .eq. 0 ) then
        if(ip_y .eq. jpe-1) then
          ien = iml-4
        else
          ien = iml-2
        endif
        if(ip_y .eq. 0) then
          isn = 4
        else
          isn = 3
        endif
        do i=isn,ien
          dummy(3)=dummy(3)+vmbd(i,4)*dx*cs(4)
        enddo
      endif
!  north
      if(ip_y .eq. jpe-1 ) then
        if(ip_y .eq. jpe-1) then
          ien = iml-4
        else
          ien = iml-2
        endif
        if(ip_y .eq. 0) then
          isn = 4
        else
          isn = 3
        endif
        do i=isn,ien
          dummy(4)=dummy(4)+vmbd(i,jml-3)*dx*cs(jml-3)
        enddo
      endif

      call sum_all(dummy,4)
      twtr=dummy(1)
      tetr=dummy(2)
      tstr=dummy(3)
      tntr=dummy(4)

! redundant
      tstr0=0.d0
      tntr0=0.d0
      twtr0=0.d0
      tetr0=0.d0
      coefx=0.5d0+0.5d0/dble(inc)
      coefy=0.5d0+0.5d0/dble(jnc)
      do j=jsm,jsm+jg-1
         twtr0=twtr0+(umom(ism-1,j)
     $        +(umom(ism,j)-umom(ism-1,j))*coefx
     $        )*dyom*texom(ism,j,1)
         tetr0=tetr0+(umom(ism+ig-1,j)
     $        +(umom(ism+ig,j)-umom(ism+ig-1,j))*coefx
     $        )*dyom*texom(ism+ig,j,1)
      enddo

      do i=ism,ism+ig-1
         tstr0=tstr0+
     $        (vmom(i,jsm-1)*csom(jsm-1)
     $        +(vmom(i,jsm)*csom(jsm)
     $        -vmom(i,jsm-1)*csom(jsm-1))*coefy
     $        )*dxom*texom(i,jsm,1)
         tntr0=tntr0+
     $        (vmom(i,jsm+jg-1)*csom(jsm+jg-1)
     $        +(vmom(i,jsm+jg)*csom(jsm+jg-1)
     $        -vmom(i,jsm+jg-1)*csom(jsm+jg))*coefy
     $        )*dxom*texom(i,jsm+jg,1)
      enddo

#ifndef ASSIM
      if(ip .eq. imaster) then
        write(*,*) 'boundary data volume daignostic'
        write(*,*) ''
        write(*,*) '(interpolated, original, difference)'
        write(*,*) 'west  ',twtr*1e-12,twtr0*1e-12,(twtr-twtr0)*1e-12
        write(*,*) 'east  ',tetr*1e-12,tetr0*1e-12,(tetr-tetr0)*1e-12
        write(*,*) 'south ',tstr*1e-12,tstr0*1e-12,(tstr-tstr0)*1e-12
        write(*,*) 'north ',tntr*1e-12,tntr0*1e-12,(tntr-tntr0)*1e-12
        write(*,*) '--------------------------------'
        write(*,*) 'net   ',(twtr-tetr+tstr-tntr)*1e-12,
     &     (twtr0-tetr0+tstr0-tntr0)*1.e-12
        write(*,*) ''
      endif
#endif
      endif


!------------------- aice ------------------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o2
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : aice'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_2d_r8(aicebd,62)
      else

#ifdef ICE
      call bcast_dble(r8o2,im0*jm0)
      call intplh(aicebd,r8o2)
c
      do j = 1, jml
      do i = 1, iml
         aicebd(i,j)=aicebd(i,j)*tex(i,j,1)
#ifdef DEBUG
         if(nint(aicebd(i,j))==999) then
            write(*,*) 'unexpected value for aicebd'
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo

#endif
      endif

      if(iflag_off .eq. 2) then
        call write_2d8_r8(aicebd,62)
      endif

!--------------------------- volice ---------------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o2
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : volice'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_2d_r8(volicebd,62)
      else

#ifdef ICE
      call bcast_dble(r8o2,im0*jm0)

!  Insted of volice, using hice for boundary conditions, do not need.
      do j=1,jm0
      do i=1,im0
         if(texom(i,j,1).ne.0.d0) then
            r8o2(i,j)=r8o2(i,j)*dxdeg*dydeg/dxdeg0/dydeg0
         endif
      enddo
      enddo

      call intplh(volicebd,r8o2)

      do j = 1, jml
      do i = 1, iml
         volicebd(i,j)=volicebd(i,j)*tex(i,j,1)
#ifdef DEBUG
         if(nint(volicebd(i,j))==999) then
            write(*,*) 'unexpected value for volicebd'
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo
#endif
      endif

      if(iflag_off .eq. 2) then
        call write_2d8_r8(volicebd,62)
      endif

#ifdef ICE
!  aice & volice limitations

      do j = 1, jml
      do i = 1, iml
      if(aicebd(i,j).lt.aicemin.or.volicebd(i,j).le.0.d0) then
         aicebd(i,j) = 0.
         volicebd(i,j) = 0.
      else
         if(aicebd(i,j).gt.1.d0) then
            aicebd(i,j) = 1.d0
         endif
         h_ice=volicebd(i,j)/(aicebd(i,j)*areat(i,j,1)+1.-tex(i,j,1))
     &           *tex(i,j,1)
         if(h_ice.lt.hmin) then
            aicebd(i,j)=volicebd(i,j)/(areat(i,j,1)*hmin+1.-tex(i,j,1))
            if(aicebd(i,j).lt.aicemin) then
               aicebd(i,j) = 0.
               volicebd(i,j) = 0.
            endif
         endif
         if(h_ice.gt.hmax) then
            aicebd(i,j)=volicebd(i,j)/(areat(i,j,1)*hmax+1.-tex(i,j,1))
            if(aicebd(i,j).gt.1.d0) then
               aicebd(i,j) = 1.d0
                  volicebd(i,j) = areat(i,j,1)*hmax
               endif
            endif
         endif
      enddo
      enddo
#endif

!---------------------------- uice ----------------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o2
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : uice'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_2d_r8(uicebd,62)
      else

#ifdef ICE
      call bcast_dble(r8o2,im0*jm0)
      call intplum(uicebd,r8o2)

      do j = 1, jml
      do i = 1, iml
         uicebd(i,j)=uicebd(i,j)*ex(i,j,1)
#ifdef DEBUG
         if(nint(uicebd(i,j))==999) then
            write(*,*) 'unexpected value for uicebd'
#ifndef DEBUG20
            stop
#endif
         endif
#endif
      enddo
      enddo
#endif
      endif

      if(iflag_off .eq. 2) then
        call write_2d8_r8(uicebd,62)
      endif

!--------------------------- vice ---------------------------------
      if(ip.eq.imaster) then
        read(omunt) r8o2
#ifdef DEBUG
        if(ionum/=0) then
           write(*,*) 'error number=',ionum
           write(*,*) 'field data reading err : vice'
           stop
        endif
#endif
      endif

      if(iflag_off .eq. 1) then
        call read_2d_r8(vicebd,62)
      else

      call bcast_dble(r8o2,im0*jm0)
#ifdef ICE
      call intplvm(vicebd,r8o2)

      do j = 1, jml
      do i = 1, iml
         vicebd(i,j)=vicebd(i,j)*ex(i,j,1)
#ifdef DEBUG
      if(nint(vicebd(i,j))==999) then
         write(*,*) 'unexpected value for vicebd'
         stop
      endif
#endif
      enddo
      enddo
#endif
      endif

      if(iflag_off .eq. 2) then
        call write_2d8_r8(vicebd,62)
      endif

!---------------- baroclinic velocities --------------------------
!  local
      do j=1,jml-1
      do i=1,iml-1
        hclubd(i,j)=ex(i,j,1)*(ashf(j)*(hclbd(i,j)+hclbd(i+1,j))
     &    +anhf(j)*(hclbd(i,j+1)+hclbd(i+1,j+1)))/areauu(j)
      enddo
      enddo

      do k=1,km
      do j=1,jml
      do i=1,iml
         dztbd(i,j,k)=dzt(i,j,k)
         dzubd(i,j,k)=dzu(i,j,k)
      enddo
      enddo
      enddo

      do k = 1,kadp
      do j = 1,jml-1
      do i = 1,iml-1
         dztbd(i,j,k)=tex(i,j,k)*(dz(k)+hclbd(i,j)*dz(k)/depadp)
         dzubd(i,j,k)=ex(i,j,k)*(dz(k)+hclubd(i,j)*dz(k)/depadp)
      enddo
      enddo
      enddo


!===== calculate baroclinic velocity with thermal wind equation ========
!------- Integrading from surface to bottom ------------------


      do k=0,km+1
      do j=1,jm
      do i=1,im
         ubd(i,j,k)=0.d0
         vbd(i,j,k)=0.d0
         rho(i,j,k)=0.d0
      enddo
      enddo
      enddo


      do k = 1,km
      do j = 1,jml
      do i = 1,iml

      if(tex(i,j,k).eq.1.d0)then

        rhoo(i,j,k)= 999.842594d0 + 6.793952d-2*tbd(i,j,k)
     &     - 9.095290d-3*tbd(i,j,k)*tbd(i,j,k)
     $     + 1.001685d-4*tbd(i,j,k)**3
     &     - 1.120083d-6*tbd(i,j,k)**4 + 6.536332d-9*tbd(i,j,k)**5
     &     + sbd(i,j,k)*(0.824493d0 - 4.0899d-3*tbd(i,j,k)
     &     + 7.6438d-5*tbd(i,j,k)*tbd(i,j,k)
     &     - 8.2467d-7*tbd(i,j,k)**3 + 5.3875d-9*tbd(i,j,k)**4)
     &     + dsqrt(sbd(i,j,k)**3)*(-5.72466d-3 + 1.0227d-4*tbd(i,j,k)
     &     - 1.6546d-6*tbd(i,j,k)*tbd(i,j,k)) + 4.8314d-4*sbd(i,j,k)**2

        rhoo(i,j,k)=rhoo(i,j,k)*1.d-3-1.d0

!      compute rho(s,theta,p)

        stbup = 1.965933d4 + 1.444304d2*tbd(i,j,k)
     $  - 1.706103d0*tbd(i,j,k)*tbd(i,j,k)
     &  + 9.648704d-3*tbd(i,j,k)**3  - 4.190253d-5*tbd(i,j,k)**4
     &  + sbd(i,j,k)*(52.84855d0 - 3.101089d-1*tbd(i,j,k)
     &  + 6.283263d-3*tbd(i,j,k)*tbd(i,j,k) -5.084188d-5*tbd(i,j,k)**3)
     &  + dsqrt(sbd(i,j,k)**3)*(3.886640d-1 + 9.085835d-3*tbd(i,j,k)
     &  - 4.619924d-4*tbd(i,j,k)*tbd(i,j,k))
     &  + pdom(k)*(3.186519d0 + 2.212276d-2*tbd(i,j,k)
     &  - 2.984642d-4*tbd(i,j,k)*tbd(i,j,k) + 1.956415d-6*tbd(i,j,k)**3)
     &  + pdom(k)*sbd(i,j,k)*(6.704388d-3  -1.847318d-4*tbd(i,j,k)
     &  + 2.059331d-7*tbd(i,j,k)*tbd(i,j,k))
     $  + 1.480266d-4*pdom(k)*dsqrt(sbd(i,j,k)**3)
     &  + pdom(k)*pdom(k)*(2.102898d-4 - 1.202016d-5*tbd(i,j,k)
     &  + 1.394680d-7*tbd(i,j,k)*tbd(i,j,k))
     $  + pdom(k)*pdom(k)*sbd(i,j,k)*(-2.040237d-6
     &  + 6.128773d-8*tbd(i,j,k) + 6.207323d-10*tbd(i,j,k)*tbd(i,j,k))

        rho(i,j,k) = ((rhoo(i,j,k)+1.d0)/(1.d0 - pdom(k)/stbup) -1.d0
     $       )*tex(i,j,k)

#ifdef DEBUG1
      if(ip.eq.0 .and. i.eq.26 .and. j.eq.1 .and. k.le.5) then
        write(*,*) 'rho,',i,j,k,rho(i,j,k),t(i,j,k),s(i,j,k),rhoo(i,j,k)
        write(*,*)
      endif
#endif

      endif
      enddo
      enddo
      enddo

      do 100 j = 1,jml-1

!     k = 1
      do i = 1,iml-1
         ubd(i,j,1)=ubd(i,j,0)
     &        -ex(i,j,1)*grav*dy2r*ddmnar*0.5d0*dzubd(i,j,1)
     $        *(rho(i+1,j+1,1)-rho(i+1,j,1)+rho(i,j+1,1)-rho(i,j,1))
     $        /cor(j)

         vbd(i,j,1)=vbd(i,j,0)
     &        +ex(i,j,1)*grav*dx2r*csr(j)*ddmnar*0.5d0*dzubd(i,j,1)
     $        *(rho(i+1,j+1,1)+rho(i+1,j,1)-rho(i,j+1,1)-rho(i,j,1))
     $        /cor(j)
#ifdef DEBUG1
      if(ip.eq.0 .and. i.eq.26 .and. j.eq.1) then
        write(*,*) 'rho,',i,j,k,dzubd(i,j,1),rho(i+1,j+1,1),
     &   rho(i+1,j,1),rho(i,j+1,1),rho(i,j,1)
        write(*,*) 'hcl,',hclbd(i+1,j+1),hclbd(i,j+1),
     &   hclbd(i,j+1),hclbd(i,j)
        write(*,*) 't,',tbd(i+1,j+1,1),tbd(i,j+1,1),
     &   tbd(i,j+1,1),tbd(i,j,1)
        write(*,*) 's,',sbd(i+1,j+1,1),sbd(i,j+1,1),
     &   sbd(i,j+1,1),sbd(i,j,1)
        write(*,*)
      endif
#endif
      enddo

      do k = 2,kadp
      do i = 1,iml-1
         ubd(i,j,k)=ubd(i,j,k-1)-ex(i,j,k)
     $        *grav*ddmnar*.5d0*(dzubd(i,j,k)+dzubd(i,j,k-1))
     $        *dy2r*( 0.5d0*(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i,j+1,k)-rho(i+1,j,k)-rho(i,j,k))
     $        +(hclbd(i+1,j+1)+hclbd(i,j+1)-hclbd(i+1,j)-hclbd(i,j))
     &        /depadp*0.25d0
     $        *(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)+rho(i+1,j,k-1)
     $        +rho(i,j,k-1)-rho(i+1,j+1,k)-rho(i,j+1,k)
     $        -rho(i+1,j,k)-rho(i,j,k)) )
     $        /cor(j)
         vbd(i,j,k)=vbd(i,j,k-1)+ex(i,j,k)
     $        *grav*ddmnar*0.5d0*(dzubd(i,j,k)+dzubd(i,j,k-1))
     $        *dx2r*( 0.5*csr(j)*(rho(i+1,j+1,k-1)+rho(i+1,j,k-1)
     $        -rho(i,j+1,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i+1,j,k)-rho(i,j+1,k)-rho(i,j,k))
     $        +(hclbd(i+1,j+1)-hclbd(i,j+1))/depadp*0.25d0
     $        *cstr(j+1)*(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j+1,k)-rho(i,j+1,k))
     $        +(hclbd(i+1,j)-hclbd(i,j))/depadp*0.25d0
     $        *cstr(j)*(rho(i+1,j,k-1)+rho(i,j,k-1)
     $        -rho(i+1,j,k)-rho(i,j,k)) )
     $        /cor(j)
#ifdef DEBUG1
      if(ip.eq.0 .and. i.eq.26 .and. j.eq.1) then
        write(*,*) 'rho,',i,j,k,dzubd(i,j,k),rho(i+1,j+1,k),
     &   rho(i+1,j,k),rho(i,j+1,k),rho(i,j,k)
        write(*,*) 't,',tbd(i+1,j+1,k),tbd(i,j+1,k),
     &   tbd(i,j+1,k),tbd(i,j,k)
        write(*,*) 's,',sbd(i+1,j+1,k),sbd(i,j+1,k),
     &   sbd(i,j+1,k),sbd(i,j,k)
        write(*,*)
      endif
#endif
      enddo
      enddo

      k = kadp+1

      do i = 1,iml-1
         ubd(i,j,k)=ubd(i,j,k-1)-ex(i,j,k)
     $        *grav*dy4r*ddmnar
     $        *(dzz(k)*.5d0*(1.d0+dzr(k)*dzubd(i,j,k))
     $        +0.5d0*hclubd(i,j)*dep(kadp)/depadp)
     $        *(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i,j+1,k)-rho(i+1,j,k)-rho(i,j,k))
     $        /cor(j)

         vbd(i,j,k)=vbd(i,j,k-1)+ex(i,j,k)*csr(j)
     $        *grav*dx4r*ddmnar
     $        *(dzz(k)*.5d0*(1.d0+dzr(k)*dzubd(i,j,k))
     $        +0.5d0*hclubd(i,j)*dep(kadp)/depadp)
     $        *(rho(i+1,j+1,k-1)+rho(i+1,j,k-1)
     $        -rho(i,j+1,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i+1,j,k)-rho(i,j+1,k)-rho(i,j,k))
     $        /cor(j)
#ifdef DEBUG1
      if(ip.eq.0 .and. i.eq.26 .and. j.eq.1) then
        write(*,*) 'rho,',i,j,k,dzubd(i,j,k),rho(i+1,j+1,k),
     &   rho(i+1,j,k),rho(i,j+1,k),rho(i,j,k)
        write(*,*) 't,',tbd(i+1,j+1,k),tbd(i,j+1,k),
     &   tbd(i,j+1,k),tbd(i,j,k)
        write(*,*) 's,',sbd(i+1,j+1,k),sbd(i,j+1,k),
     &   sbd(i,j+1,k),sbd(i,j,k)
        write(*,*)
      endif
#endif

      enddo

      do k = kadp+2,km
      do i = 1,iml-1
         ubd(i,j,k)=ubd(i,j,k-1)-ex(i,j,k)
     $        *grav*dy4r*dzz(k)*ddmnar*.5d0*(1.d0+dzr(k)*dzubd(i,j,k))
     $        *(rho(i+1,j+1,k-1)+rho(i,j+1,k-1)
     $        -rho(i+1,j,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i,j+1,k)-rho(i+1,j,k)-rho(i,j,k))
     $        /cor(j)

         vbd(i,j,k)=vbd(i,j,k-1)+ex(i,j,k)*csr(j)
     $        *grav*dx4r*dzz(k)*ddmnar*.5d0*(1.d0+dzr(k)*dzubd(i,j,k))
     $        *(rho(i+1,j+1,k-1)+rho(i+1,j,k-1)
     $        -rho(i,j+1,k-1)-rho(i,j,k-1)
     $        +rho(i+1,j+1,k)+rho(i+1,j,k)-rho(i,j+1,k)-rho(i,j,k))
     $        /cor(j)
#ifdef DEBUG1
      if(ip.eq.0 .and. i.eq.26 .and. j.eq.1) then
        write(*,*) 'rho,',i,j,k,dzubd(i,j,k),rho(i+1,j+1,k),
     &   rho(i+1,j,k),rho(i,j+1,k),rho(i,j,k)
        write(*,*) 't,',tbd(i+1,j+1,k),tbd(i,j+1,k),
     &   tbd(i,j+1,k),tbd(i,j,k)
        write(*,*) 's,',sbd(i+1,j+1,k),sbd(i,j+1,k),
     &   sbd(i,j+1,k),sbd(i,j,k)
        write(*,*)
      endif
#endif

      enddo
      enddo

      do i=1,iml-1
         bttu=0.d0
         bttv=0.d0
         do k=1,km
            bttu=bttu + ubd(i,j,k)*dzubd(i,j,k)
     &           /(hrr(i,j)+hclubd(i,j)+(1.-ex(i,j,k)))*ex(i,j,k)
            bttv=bttv + vbd(i,j,k)*dzubd(i,j,k)
     &           /(hrr(i,j)+hclubd(i,j)+(1.-ex(i,j,k)))*ex(i,j,k)
         enddo

         do k=1,km
            ubd(i,j,k)=(ubd(i,j,k)-bttu)*ex(i,j,k)
            vbd(i,j,k)=(vbd(i,j,k)-bttv)*ex(i,j,k)
#ifdef DEBUG1
      if(ip.eq.0 .and. i.eq.26 .and. j.eq.1) then
        write(*,*) i,j,k,ubd(i,j,k),bttu,ubd(i,j,k)+bttu
        write(*,*) i,j,k,vbd(i,j,k),bttv,vbd(i,j,k)+bttv
        write(*,*) dzubd(i,j,k),hrr(i,j),hclubd(i,j)
        write(*,*)
      endif
#endif
         enddo
      enddo

 100  continue


      do j=1,jml-1
      do i=1,iml-1
         sfunbd(i,j)=ex(i,j,1)*umbd(i,j)
     &        /(hrr(i,j)+hclubd(i,j)+(1.d0-ex(i,j,1)))
         sfvnbd(i,j)=ex(i,j,1)*vmbd(i,j)
     &        /(hrr(i,j)+hclubd(i,j)+(1.d0-ex(i,j,1)))
      enddo
      enddo

#ifdef NWNPAC
      do k=1,km
      do j=1,jml
      do i=1,iml
         tref12(i,j,k,m)=tbd(i,j,k)
         sref12(i,j,k,m)=sbd(i,j,k)
      enddo
      enddo
      enddo
#endif
#ifdef JP68M
      do k=1,km
      do j=1,jml
      do i=1,iml
         tref12(i,j,k,m)=tbd(i,j,k)
         sref12(i,j,k,m)=sbd(i,j,k)
      enddo
      enddo
      enddo
#endif
#ifdef ROKKA
      do k=1,km
      do j=1,jml
      do i=1,iml
         tref12(i,j,k,m)=tbd(i,j,k)
         sref12(i,j,k,m)=sbd(i,j,k)
      enddo
      enddo
      enddo
#endif

      if (m.eq.1) then
!----------------- initial conditions ----------------------

      if (nfirst) then

#ifndef JP44

#if defined(ROKKA) && defined(EXP2006)
      if(ip==imaster .and. num_inidt>1) then
         do l=1,14*(num_inidt-1)
            read(contini)
         enddo
      endif
#endif

      if(ip.eq.imaster) then
         read(contini) idum,adum
         read(contini) pd,pm,ddmna,dmn
      endif
!      call bcast_int1(idum,1)
      call bcast_dble1(adum,1)
!      ctime_ini=adum*3.6d3*dble(nday)/dble(nday_obc)
!      ctime=ctime_ini-dtts
!      nsec=adum*dble(nday)/dble(nday_obc)*36.d2

      call bcast_dble(pd,km+2)
      call bcast_dble(pm,km+2)
      call bcast_dble1(ddmna,1)
      call bcast_dble(dmn,km+2)

      adum=adum*dble(nday)/dble(nday_obc)*36.d2
!      if(obc_clim.ne.1) ahour=adum
#if defined(EXP2003)
      adum=adum-dble(360*2*86400)
#elif defined(EXP2006) || defined(EXP2007)
      do while(adum>nday*86400)
         adum=adum-dble(nday*86400)
      enddo
#endif
      nsec=int(adum)+1
#ifdef ROKKA
      nsec=nsec+inihour*3600
#endif
      if(nsec>=nday*86400) nsec=nsec-nday*86400

      if(ip==imaster)write(*,*)'nsec=',nsec


! ---- read u ----

      if(ip.eq.imaster) then
#ifdef ROKKA
         read(contini) r8o3
#else
!#ifdef CASE51
         read(contini) r8o3
!#else
!         do k=1,km
!            read(contini) ((r8o3(i,j,k),i=1,im0),j=1,jm0)
!         enddo
!#endif
#endif
      endif

#ifndef REST_INI
      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplu(u,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
         u(i,j,k)=u(i,j,k)*ex(i,j,k)
#ifdef DEBUG
         if(nint(u(i,j,k))==999) then
            write(*,*) 'unexpected value for initial u at'
            write(*,*) '  i=',i-1+lw(ip_x),', j=',j-1+ls(ip_y),', k=',k
            stop
         endif
#endif
      enddo
      enddo
      enddo
#else
      do k=1,km
      do j=1,jml
      do i=1,iml
         u(i,j,k)=0.d0
      enddo
      enddo
      enddo
#endif
! ---- read v ----

      if(ip.eq.imaster) then
#ifdef ROKKA
         read(contini) r8o3
#else
!#ifdef CASE51
         read(contini) r8o3
!#else
!         do k=1,km
!            read(contini) ((r8o3(i,j,k),i=1,im0),j=1,jm0)
!         enddo
!#endif
#endif
      endif


#ifndef REST_INI
      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplu(v,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
         v(i,j,k)=v(i,j,k)*ex(i,j,k)
#ifdef DEBUG
         if(nint(v(i,j,k))==999) then
            write(*,*) 'unexpected value for initial v at'
            write(*,*) '  i=',i-1+lw(ip_x),', j=',j-1+ls(ip_y),', k=',k
            stop
         endif
#endif
      enddo
      enddo
      enddo
#else
      do k=1,km
      do j=1,jml
      do i=1,iml
         v(i,j,k)=0.d0
      enddo
      enddo
      enddo
#endif


!---- read t ----

      if(ip.eq.imaster) then
#ifdef ROKKA
         read(contini) r8o3
#else
!#ifdef CASE51
         read(contini) r8o3
!#else
!         do k=1,km
!            read(contini) ((r8o3(i,j,k),i=1,im0),j=1,jm0)
!         enddo
!#endif
#endif
      endif

      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplt(t,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
         t(i,j,k)=t(i,j,k)*tex(i,j,k)
#ifdef DEBUG
         if(nint(t(i,j,k))==999) then
            write(*,*) 'unexpected value for initial t at'
            write(*,*) '  i=',i-1+lw(ip_x),', j=',j-1+ls(ip_y),', k=',k
            stop
         endif
#endif
      enddo
      enddo
      enddo

!---- read s ----

      if(ip.eq.imaster) then
#ifdef ROKKA
         read(contini) r8o3
#else
!#ifdef CASE51
         read(contini) r8o3
!#else
!         do k=1,km
!            read(contini) ((r8o3(i,j,k),i=1,im0),j=1,jm0)
!         enddo
!#endif
#endif
      endif

      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplt(s,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
         s(i,j,k)=s(i,j,k)*tex(i,j,k)
#ifdef DEBUG
         if(nint(s(i,j,k))==999) then
            write(*,*) 'unexpected value for initial s at'
            write(*,*) '  i=',i-1+lw(ip_x),', j=',j-1+ls(ip_y),', k=',k
            stop
         endif
#endif
      enddo
      enddo
      enddo

!---- read tke ----

      if(ip.eq.imaster) then
#ifdef ROKKA
         read(contini) r8o3
#else
!#ifdef CASE51
         read(contini) r8o3
!#else
!    do k=1,km
!            read(contini) ((r8o3(i,j,k),i=1,im0),j=1,jm0)
!         enddo
!#endif
#endif
      endif

      call bcast_dble(r8o3,im0*jm0*(km+2))
      call intplt(tke,r8o3)

      do k = 1, km
      do j = 1, jml
      do i = 1, iml
         tke(i,j,k)=tke(i,j,k)*tex(i,j,k)
#ifdef DEBUG
         if(tke(i,j,k)==999.d0) then
            write(*,*) 'unexpected value for initial s at'
            write(*,*) '  i=',i-1+lw(ip_x),', j=',j-1+ls(ip_y),', k=',k
            stop
         endif
#endif
      enddo
      enddo
      enddo

!---- read hcl ----

      if(ip.eq.imaster) then
        read(contini) r8o2
      endif

      call bcast_dble(r8o2,im0*jm0)
      call intplh(hcl,r8o2)

      do j = 1, jml
      do i = 1, iml
         hcl(i,j)=hcl(i,j)*tex(i,j,1)
#ifdef DEBUG
         if(hcl(i,j)==999.d0) then
            write(*,*) 'unexpected value for hcl'
            stop
         endif
#endif
      enddo
      enddo
!#endif
!---- read um ----

      if(ip.eq.imaster) then
        read(contini) r8o2
      endif

#ifndef REST_INI
      call bcast_dble(r8o2,im0*jm0)
      call intplum(um,r8o2)

      do j = 1, jml
      do i = 1, iml
         um(i,j)=um(i,j)*ex(i,j,1)
#ifdef DEBUG
      if(um(i,j)==999.d0) then
         write(*,*) 'unexpected value for initial um'
         stop
      endif
#endif
      enddo
      enddo

#else
      do j=1,jml
      do i=1,iml
         um(i,j)=0.d0
      enddo
      enddo
#endif

!---- read vm ----

      if(ip.eq.imaster) then
        read(contini) r8o2
      endif

#ifndef REST_INI
      call bcast_dble(r8o2,im0*jm0)
      call intplvm(vm,r8o2)

      do j = 1, jml
      do i = 1, iml
         vm(i,j)=vm(i,j)*ex(i,j,1)
#ifdef DEBUG
      if(vm(i,j)==999.d0) then
         write(*,*) 'unexpected value for initial vm'
         stop
      endif
#endif
      enddo
      enddo
#else
      do j=1,jml
      do i=1,iml
         vm(i,j)=0.d0
      enddo
      enddo
#endif

!---- read aice ----

      if(ip.eq.imaster) then
        read(contini) r8o2
      endif

#ifndef ICE
      do j=1,jml
      do i=1,iml
         aice(i,j)=0.d0
      enddo
      enddo
#else
      call bcast_dble(r8o2,im0*jm0)
      call intplh(aice,r8o2)

      do j = 1, jml
      do i = 1, iml
         aice(i,j)=aice(i,j)*tex(i,j,1)
#ifdef DEBUG
         if(nint(aice(i,j))==999) then
            write(*,*) 'unexpected value for initial aice'
            stop
         endif
#endif
      enddo
      enddo
#endif

!---- read volice ----

      if(ip.eq.imaster) then
        read(contini) r8o2
      endif

#ifndef ICE
      do j=1,jml
      do i=1,iml
         volice(i,j)=0.d0
      enddo
      enddo
#else
      call bcast_dble(r8o2,im0*jm0)
      call intplh(volice,r8o2)

      do j = 1, jml
      do i = 1, iml
#ifdef DEBUG
         if(texn(i,j)>0 .and.nint(volice(i,j))==999) then
            write(*,*) 'unexpected value for initial volice'
            stop
         endif
#endif
         volice(i,j)=volice(i,j)*dxdeg*dydeg/dxdeg0/dydeg0*tex(i,j,1)
      enddo
      enddo
#endif

!---- read uice ----

      if(ip.eq.imaster) then
        read(contini) r8o2
      endif

#if !defined(ICE) || defined(REST_INI)
      do j=1,jml
      do i=1,iml
         uice(i,j)=0.d0
      enddo
      enddo

#else
      call bcast_dble(r8o2,im0*jm0)
      call intplum(uice,r8o2)

      do j = 1, jml
      do i = 1, iml
         uice(i,j)=uice(i,j)*ex(i,j,1)
#ifdef DEBUG
         if(nint(uice(i,j))==999) then
            write(*,*) 'unexpected value for initial uice'
            stop
         endif
#endif
      enddo
      enddo
#endif

!---- read vice ----

      if(ip.eq.imaster) then
        read(contini) r8o2
      endif

#if !defined(ICE) || defined(REST_INI)
      do j=1,jml
      do i=1,iml
         vice(i,j)=0.d0
      enddo
      enddo
#else
      call bcast_dble(r8o2,im0*jm0)
      call intplvm(vice,r8o2)

      do j = 1, jml
      do i = 1, iml
         vice(i,j)=vice(i,j)*ex(i,j,1)
#ifdef DEBUG
         if(nint(vice(i,j))==999) then
            write(*,*) 'unexpected value for initial vice'
            stop
         endif
#endif
      enddo
      enddo
#endif


#else

      nsec=int(ahourbd(1)*3600.)+1
      do while(nsec>=nday*86400)
         nsec=nsec-nday*86400
      enddo
      if(nsec>=nday*86400) nsec=nsec-nday*86400


      do k=1,km
      do j=1,jml
      do i=1,iml
         t(i,j,k)=tbd(i,j,k)*tex(i,j,k)
         s(i,j,k)=sbd(i,j,k)*tex(i,j,k)
         tke(i,j,k)=tkebd(i,j,k)*tex(i,j,k)
#ifdef REST_INI
         u(i,j,k)=0.d0
         v(i,j,k)=0.d0
#else
         u(i,j,k)=(ubd(i,j,k)+sfunbd(i,j))*ex(i,j,k)
         v(i,j,k)=(vbd(i,j,k)+sfvnbd(i,j))*ex(i,j,k)
#endif
      enddo
      enddo
      enddo

      do j=1,jml
      do i=1,iml
         hcl(i,j)=hclbd(i,j)*tex(i,j,1)
         hclu(i,j)=hclubd(i,j)*ex(i,j,1)
#ifdef REST_INI
         um(i,j)=0.d0
         vm(i,j)=0.d0
         uice(i,j)=0.d0
         vice(i,j)=0.d0
#else
         um(i,j)=umbd(i,j)*ex(i,j,1)
         vm(i,j)=vmbd(i,j)*ex(i,j,1)
         uice(i,j)=uicebd(i,j)*ex(i,j,1)
         vice(i,j)=vicebd(i,j)*ex(i,j,1)
#endif
         aice(i,j)=aicebd(i,j)*tex(i,j,1)
         volice(i,j)=volicebd(i,j)*tex(i,j,1)
      enddo
      enddo
#endif
      else
      if(ip.eq.imaster) then
!----------------- temporaly -----------------------
!        read(restart) last,ahour
!---------------------------------------------------
         read(restart) last,ahour,nsec
      endif

#ifdef ASSIM
! exp0201 start from 02/01/13
!        ahourini=24.*dble(12+30*0+360*2)
! exp0202 start from 02/02/17
!        ahourini=24.*dble(17+30*1+360*2)
! exp0203 start from 02/03/24
!        ahourini=24.*dble(24+30*2+360*2)
! exp0204 start from 02/04/28
!        ahourini=24.*dble(28+30*3+360*2)
! exp0205 start from 02/06/02
!        ahourini=24.*dble(2+30*5+360*2)
! exp0206 start from 02/07/07
!        ahourini=24.*dble(7+30*6+360*2)
! exp0207 start from 02/08/11
!        ahourini=24.*dble(11+30*7+360*2)
! exp0208 start from 02/09/15
!        ahourini=24.*dble(15+30*8+360*2)
! exp020 start from 02/010/20
!        ahourini=24.*dble(20+30*9+360*2)
! exp0210 start from 02/11/24
!        ahourini=24.*dble(24+30*10+360*2)
! exp0301 start from 02/02/07
!        ahourini=24.*dble(7+30*1+360*2)
! exp0302 start from 02/03/07
!        ahourini=24.*dble(7+30*2+360*2)
! exp0303 start from 02/04/04
!        ahourini=24.*dble(4+30*3+360*2)
! exp0304 start from 02/05/02
!        ahourini=24.*dble(2+30*4+360*2)
! exp0305 start from 02/05/30
!        ahourini=24.*dble(30+30*4+360*2)
! exp0306 start from 02/06/27
!        ahourini=24.*dble(27+30*5+360*2)
! exp0307 start from 02/07/25
!        ahourini=24.*dble(25+30*6+360*2)
! exp0308 start from 02/08/22
!        ahourini=24.*dble(22+30*7+360*2)
! exp0309 start from 02/09/19
!        ahourini=24.*dble(19+30*8+360*2)
! exp0310 start from 02/10/17
!        ahourini=24.*dble(17+30*9+360*2)
! exp0311 start from 02/11/14
!        ahourini=24.*dble(14+30*10+360*2)
! exp0312 start from 02/12/12
!        ahourini=24.*dble(12+30*11+360*2)
! exp0313 start from 02/01/09
!        ahourini=24.*dble(9+30*0+360*2)
! exp0314 start from 02/02/06
!        ahourini=24.*dble(6+30*1+360*2)
! exp0315 start from 02/03/05
!        ahourini=24.*dble(5+30*2+360*2)
! exp0316 start from 02/04/02
!        ahourini=24.*dble(2+30*3+360*2)
! exp0317 start from 02/04/30
!        ahourini=24.*dble(30+30*3+360*2)
! exp0318 start from 02/05/28
!        ahourini=24.*dble(28+30*4+360*2)

! exp0503 start from 02/04/08
        ahourini=24.*dble(8+30*3+360*2)
! exp0504 start from 02/05/06
!        ahourini=24.*dble(6+30*4+360*2)
! exp0505 start from 02/06/3
!        ahourini=24.*dble(3+30*5+360*2)
! exp0506 start from 02/07/1
!        ahourini=24.*dble(1+30*6+360*2)
! exp0507 start from 02/07/29
!        ahourini=24.*dble(29+30*6+360*2)
! exp0508 start from 02/08/26
!        ahourini=24.*dble(26+30*7+360*2)
! exp0509 start from 02/09/23
!        ahourini=24.*dble(23+30*8+360*2)
! exp0510 start from 02/10/20
!        ahourini=24.*dble(20+30*9+360*2)


        nkaiini=nint(ahourini*60./dtuv)
! ctime is always between year 80 and year 81
!        ctini=(ahourini+360.d0*78.d0*24.d0)*3600.
        nsecini= mod( nint(ahourini*3600.),nday*86400)
        nsec = nsecini
!        ctime=ctini-dtts
#else
      call bcast_int1(last,1)
      call bcast_dble1(ahour,1)
!-------------- temporaly -----------------
!      nsec=1
!------------------------------------------
      call bcast_int1(nsec,1)
!------------------------------------------
!      call bcast_dble1(ctime_ini,1)
!      call bcast_dble1(ctime,1)
#endif
!      ctime_ini=ahourbd(1)*3.6d3 ! temporaly
!      ctime=ctime_ini+ahour*3.6d3-dtts

!      ctime=ctime-dtts

#ifndef CLIMAT
      if(nday==365 .and. nsec>nday*86400) nsec=nsec-86400
      call calender(cal_year,nsec,month,iday,ihour,imin,isec)
      if(ip==imaster) then
         write(*,*) 'start from'
         write(*,'("date ",i4,i2.2,i2.2)') cal_year,month,iday
         write(*,'("time ",i2.2,":",i2.2," ",i2.2)') ihour,imin,isec
         write(*,*)
      endif
#endif

      if(ip.eq.imaster) then
        read(restart) pd,pm,ddmna,dmn
      endif
      call bcast_dble(pd,km+2)
      call bcast_dble(pm,km+2)
      call bcast_dble1(ddmna,1)
      call bcast_dble(dmn,km+2)

      call read_3d_r8(u,restart)
      call read_3d_r8(v,restart)
      call read_3d_r8(t,restart)
      call read_3d_r8(s,restart)
      call read_3d_r8(tke,restart)

      call read_2d_r8(hcl,restart)
      call read_2d_r8(um,restart)
      call read_2d_r8(vm,restart)

      call read_2d_r8(aice,restart)
      call read_2d_r8(volice,restart)
      call read_2d_r8(uice,restart)
      call read_2d_r8(vice,restart)

      do j=1,jml-1
      do i=1,iml-1
        hclu(i,j)=ex(i,j,1)*(ashf(j)*(hcl(i,j)+hcl(i+1,j))
     &    +anhf(j)*(hcl(i,j+1)+hcl(i+1,j+1)))/areauu(j)
      enddo
      enddo

      endif

      endif

!------------------ size boundary condition -----------------------

!west
      if( ip_x .eq. 0) then
      do k = 1,km
      do j = 1,jml
         do i=1,indg+2
            twbc(i,j,k,m)=tbd(i+2,j,k)
            swbc(i,j,k,m)=sbd(i+2,j,k)
            tkewbc(i,j,k,m)=tkebd(i+2,j,k)
         enddo
         do i=1,2
            uwbc(i,j,k,m)=ubd(i+2,j,k)
            vwbc(i,j,k,m)=vbd(i+2,j,k)
         enddo
      enddo
      enddo
      do j=1,jml
        do i=1,indg+2
           hclwbc(i,j,m)=hclbd(i+2,j)
           wtwbc(i,j,m)=wtbd(i+2,j)
           volicewbc(i,j,m)=volicebd(i+2,j)
        enddo
         do i=1,2
            umwbc(i,j,m)=umbd(i+2,j)
            vmwbc(i,j,m)=vmbd(i+2,j)
            uicewbc(i,j,m)=uicebd(i+2,j)
            vicewbc(i,j,m)=vicebd(i+2,j)
            aicewbc(i,j,m)=aicebd(i+2,j)
         enddo
      enddo
      endif
!east
      if( ip_x .eq. ipe-1) then
      do k=1,km
      do j=1,jml
         do i=1,indg+2
            tebc(i,j,k,m)=tbd(iml-1-i,j,k)
            sebc(i,j,k,m)=sbd(iml-1-i,j,k)
            tkeebc(i,j,k,m)=tkebd(iml-1-i,j,k)
         enddo
         do i=1,2
            uebc(i,j,k,m)=ubd(iml-2-i,j,k)
            vebc(i,j,k,m)=vbd(iml-2-i,j,k)
         enddo
      enddo
      enddo
      do j=1,jml
         do i=1,indg+2
           hclebc(i,j,m)=hclbd(iml-1-i,j)
           wtebc(i,j,m)=wtbd(iml-1-i,j)
           voliceebc(i,j,m)=volicebd(iml-1-i,j)
        enddo
         do i=1,2
            umebc(i,j,m)=umbd(iml-2-i,j)
            vmebc(i,j,m)=vmbd(iml-2-i,j)
            uiceebc(i,j,m)=uicebd(iml-2-i,j)
            viceebc(i,j,m)=vicebd(iml-2-i,j)
            aiceebc(i,j,m)=aicebd(iml-2-i,j)
         enddo
      enddo
      endif

!  south
      if( ip_y .eq. 0) then
      do k = 1,km
      do i = 1,iml
         do j=1,jndg+2
            tsbc(i,j,k,m)=tbd(i,j+2,k)
            ssbc(i,j,k,m)=sbd(i,j+2,k)
            tkesbc(i,j,k,m)=tkebd(i,j+2,k)
         enddo
         do j=1,2
            usbc(i,j,k,m)=ubd(i,j+2,k)
            vsbc(i,j,k,m)=vbd(i,j+2,k)
         enddo
      enddo
      enddo
      do i = 1,iml
         do j=1,jndg+2
            hclsbc(i,j,m)=hclbd(i,j+2)
            wtsbc(i,j,m)=wtbd(i,j+2)
            volicesbc(i,j,m)=volicebd(i,j+2)
         enddo
         do j=1,2
            umsbc(i,j,m)=umbd(i,j+2)
            vmsbc(i,j,m)=vmbd(i,j+2)
            uicesbc(i,j,m)=uicebd(i,j+2)
            vicesbc(i,j,m)=vicebd(i,j+2)
            aicesbc(i,j,m)=aicebd(i,j+2)
         enddo
      enddo
      endif
!  north
      if(ip_y .eq. jpe-1) then
      do k=1,km
      do i=1,iml
         do j=1,jndg+2
            tnbc(i,j,k,m)=tbd(i,jml-1-j,k)
            snbc(i,j,k,m)=sbd(i,jml-1-j,k)
            tkenbc(i,j,k,m)=tkebd(i,jml-1-j,k)
         enddo
         do j=1,2
            unbc(i,j,k,m)=ubd(i,jml-2-j,k)
            vnbc(i,j,k,m)=vbd(i,jml-2-j,k)
         enddo
      enddo
      enddo

      do i=1,iml
          do j=1,jndg+2
            hclnbc(i,j,m)=hclbd(i,jml-1-j)
            wtnbc(i,j,m)=wtbd(i,jml-1-j)
            volicenbc(i,j,m)=volicebd(i,jml-1-j)
         enddo
         do j=1,2
            umnbc(i,j,m)=umbd(i,jml-2-j)
            vmnbc(i,j,m)=vmbd(i,jml-2-j)
            uicenbc(i,j,m)=uicebd(i,jml-2-j)
            vicenbc(i,j,m)=vicebd(i,jml-2-j)
            aicenbc(i,j,m)=aicebd(i,jml-1-j)
         enddo
      enddo
      endif

!------------------------------------------------------------


 300  continue

 1    format('error count =',i7)

      call sbnum

! west
      if(ip_x .eq. 0) then
      do k=1,km
      do j=1,jml
        do i=1,indg+2
         twbt(i,j,k)=(twbc(i,j,k,mb1)
     $        +(twbc(i,j,k,mb2)-twbc(i,j,k,mb1))*cbf)*tex(i+2,j,k)
         swbt(i,j,k)=(swbc(i,j,k,mb1)
     $        +(swbc(i,j,k,mb2)-swbc(i,j,k,mb1))*cbf)*tex(i+2,j,k)
         tkewbt(i,j,k) = (tkewbc(i,j,k,mb1)
     $        + (tkewbc(i,j,k,mb2)-tkewbc(i,j,k,mb1)) * cbf
     $        )*tex(i+2,j,k)
        enddo

      enddo
      enddo
      do j = 1,jml
        do i=1,indg+2
         hclwbt(i,j)=(hclwbc(i,j,mb1)
     $        +(hclwbc(i,j,mb2)-hclwbc(i,j,mb1))*cbf)*tex(i+2,j,1)
         hclbwbt(i,j)=hclwbt(i,j)
         wtwbt(i,j)=(wtwbc(i,j,mb1)
     $        +(wtwbc(i,j,mb2)-wtwbc(i,j,mb1))*cbf)*tex(i+2,j,1)
         wtbwbt(i,j)=wtwbt(i,j)
         volicewbt(i,j)=volicewbc(i,j,mb2)
     $        +(volicewbc(i,j,mb2)-volicewbc(i,j,mb1))*cbf
        enddo
         do i=1,2
         uicewbt(i,j)=(uicewbc(i,j,mb1)
     $           +(uicewbc(i,j,mb2)-uicewbc(i,j,mb1))*cbf)*ex(i+2,j,1)
         vicewbt(i,j)=(vicewbc(i,j,mb1)
     $           +(vicewbc(i,j,mb2)-vicewbc(i,j,mb1))*cbf)*ex(i+2,j,1)
         enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do k=1,km
      do j=1,jml
      do i=1,indg+2
         tebt(i,j,k)=(tebc(i,j,k,mb1)
     $        +(tebc(i,j,k,mb2)-tebc(i,j,k,mb1))*cbf)*tex(iml-1-i,j,k)
         sebt(i,j,k) = (sebc(i,j,k,mb1)
     $        +(sebc(i,j,k,mb2)-sebc(i,j,k,mb1))*cbf)*tex(iml-1-i,j,k)
         tkeebt(i,j,k) = (tkeebc(i,j,k,mb1)
     $        + (tkeebc(i,j,k,mb2)-tkeebc(i,j,k,mb1)) * cbf
     $        ) * tex(iml-1-i,j,k)
      enddo
      enddo
      enddo
      do j=1,jml
      do i=1,indg+2
         hclebt(i,j)=(hclebc(i,j,mb1)
     $        +(hclebc(i,j,mb2)-hclebc(i,j,mb1))*cbf)*tex(iml-1-i,j,1)
         hclbebt(i,j)=hclebt(i,j)
         wtebt(i,j)=(wtebc(i,j,mb1)
     $        +(wtebc(i,j,mb2)-wtebc(i,j,mb1))*cbf)*tex(iml-1-i,j,1)
         wtbebt(i,j)=wtebt(i,j)
         voliceebt(i,j)=voliceebc(i,j,mb2)
     $        +(voliceebc(i,j,mb2)-voliceebc(i,j,mb1))*cbf
        enddo
         do i=1,2
         uiceebt(i,j)=(uiceebc(i,j,mb1)
     $     +(uiceebc(i,j,mb2)-uiceebc(i,j,mb1))*cbf)*ex(iml-2-i,j,1)
         viceebt(i,j)=(viceebc(i,j,mb1)
     $     +(viceebc(i,j,mb2)-viceebc(i,j,mb1))*cbf)*ex(iml-2-i,j,1)
         enddo
      enddo
      endif
!  south
      if(ip_y .eq. 0) then
      do i=1,im
      do k=1,km
      do j=1,jndg+2
         tsbt(i,j,k)=(tsbc(i,j,k,mb1)
     $        +(tsbc(i,j,k,mb2)-tsbc(i,j,k,mb1))*cbf)*tex(i,j+2,k)
         ssbt(i,j,k)=(ssbc(i,j,k,mb1)
     $        +(ssbc(i,j,k,mb2)-ssbc(i,j,k,mb1))*cbf)*tex(i,j+2,k)
         tkesbt(i,j,k) = (tkesbc(i,j,k,mb1)
     $        + (tkesbc(i,j,k,mb2)-tkesbc(i,j,k,mb1)) * cbf
     $        )*tex(i,j+2,k)
      enddo
      enddo
      enddo
      do i=1,iml
      do j=1,jndg+2
         hclsbt(i,j)=(hclsbc(i,j,mb1)
     $        +(hclsbc(i,j,mb2)-hclsbc(i,j,mb1))*cbf)*tex(i,j+2,1)
         hclbsbt(i,j)=hclsbt(i,j)
         wtsbt(i,j)=(wtsbc(i,j,mb1)
     $           +(wtsbc(i,j,mb2)-wtsbc(i,j,mb1))*cbf)*tex(i,j+2,1)
         wtbsbt(i,j)=wtsbt(i,j)
         volicesbt(i,j)=volicesbc(i,j,mb2)
     $        +(volicesbc(i,j,mb2)-volicesbc(i,j,mb1))*cbf
        enddo
         do j=1,2
         uicesbt(i,j)=(uicesbc(i,j,mb1)
     $         +(uicesbc(i,j,mb2)-uicesbc(i,j,mb1))*cbf)*ex(i,j+2,1)
         vicesbt(i,j)=(vicesbc(i,j,mb1)
     $         +(vicesbc(i,j,mb2)-vicesbc(i,j,mb1))*cbf)*ex(i,j+2,1)
         enddo
      enddo
      endif
!  north
      if(ip_y .eq. jpe -1) then
      do k=1,km
      do i=1,iml
      do j=1,jndg+2
         tnbt(i,j,k)=(tnbc(i,j,k,mb1)
     $        +(tnbc(i,j,k,mb2)-tnbc(i,j,k,mb1))*cbf)*tex(i,jml-1-j,k)
         snbt(i,j,k)=(snbc(i,j,k,mb1)
     $        +(snbc(i,j,k,mb2)-snbc(i,j,k,mb1))*cbf)*tex(i,jml-1-j,k)
         tkenbt(i,j,k) = (tkenbc(i,j,k,mb1)
     $        + (tkenbc(i,j,k,mb2)-tkenbc(i,j,k,mb1)) * cbf
     $        )*tex(i,jml-1-j,k)

      enddo
      enddo
      enddo

      do i=1,iml
      do j=1,jndg+2
         hclnbt(i,j)=(hclnbc(i,j,mb1)
     $        +(hclnbc(i,j,mb2)-hclnbc(i,j,mb1))*cbf)*tex(i,jml-1-j,1)
         hclbnbt(i,j)=hclnbt(i,j)
         wtnbt(i,j)=(wtnbc(i,j,mb1)
     $        +(wtnbc(i,j,mb2)-wtnbc(i,j,mb1))*cbf)*tex(i,jml-1-j,1)
         wtbnbt(i,j)=wtnbt(i,j)
         volicenbt(i,j)=volicenbc(i,j,mb2)
     $        +(volicenbc(i,j,mb2)-volicenbc(i,j,mb1))*cbf
         enddo
         do j=1,2
         uicenbt(i,j)=(uicenbc(i,j,mb1)
     $       +(uicenbc(i,j,mb2)-uicenbc(i,j,mb1))*cbf)*ex(i,jml-2-j,1)
         vicenbt(i,j)=(vicenbc(i,j,mb1)
     $       +(vicenbc(i,j,mb2)-vicenbc(i,j,mb1))*cbf)*ex(i,jml-2-j,1)
         enddo
      enddo
      endif

#endif
      return
      end

#ifdef NESTED
!=============== interpolation for t, s and tke ==========================

      subroutine intplt(pa1,pa0)

      use param
      use mod_mpi
      implicit none
#include "common.h"

      real(8),intent(in) :: pa0(im0,jm0,0:km+1)
      real(8),intent(out) :: pa1(im,jm,0:km+1)
      real(8) :: x1,y1,big,sml
!      double precision pa1,x1,y1,big,sml
!      real*4 pa0
!      integer,parameter :: nibuf3 = 10,nibuf4=80
      real(8) ::  sp2(im,jm)
      integer :: i,j,k,nx,ny,isum,n,ii,jj

      if(inc.eq.1 .and. jnc.eq.1) then
         do k=1,km
         do j=1,jml
         do i=1,iml
            ii=lw(ip_x)+i+ism-5
            jj=ls(ip_y)+j+jsm-5
            pa1(i,j,k)=pa0(ii,jj,k)*texom(ii,jj,k)
     $           +(1.d0-texom(ii,jj,k))*999.d0

         enddo
         enddo
         enddo

         return
      endif

      big=0.9d35
      sml=-0.9d35

      do 10 k = 1, km

        nx=im
        ny=jm

         isum=isumt(k)

         do n=1,isum
            xp(n) = slon0+(ipt(n,k)-1)*dxdeg0
            yp(n) = slat0+(jpt(n,k)-1)*dydeg0
            zp(n) = pa0(ipt(n,k),jpt(n,k),k)
         enddo

         do n = isum+1, ndpt
            xp(n) = xp(isum)
            yp(n) = yp(isum)
            zp(n) = zp(isum)
         enddo

         if(isum.ne.0) then

            x1=slon+dxdeg*dble(lw(ip_x)-1)
            y1=slat+dydeg*dble(ls(ip_y)-1)

#ifdef DEBUG1
      if(ip .eq. 0) then
      write(*,*) 'level',k,' start'
      endif
#endif
#ifdef DEBUG5
      if(ip .eq. 0) then
      write(*,*) 'level',k,' start'
      endif
#endif
            call zgrid(sp2,xp,yp,zp,nx,ny,x1,y1,dxdeg,dydeg,
     $           cay,nrng,ndpt,0)
            call smooth(sp2,nx,ny,nsm)

            do j=1,ny
            do i=1,nx
               if(sp2(i,j).gt.big.or.sp2(i,j).lt.sml) sp2(i,j)=999.d0
            enddo
            enddo

            do j=1,jml
            do i=1,iml
                pa1(i,j,k)=sp2(i,j)
            enddo
            enddo
         else

            do j = 1, jml
            do i = 1, iml
               pa1(i,j,k) = 999.d0
            enddo
            enddo

         endif

 10   continue

      return
      end subroutine intplt



!================= interpolation for hcl ===============================

      subroutine intplh(pa1,pa0)

      use param
      use mod_mpi
      implicit none
#include "common.h"

      real(8),intent(in) ::  pa0(im0,jm0)
      real(8),intent(out) :: pa1(im,jm)
      real(8) ::  sp2(im,jm),x1,y1,big,sml
      integer :: i,j,n,nx,ny,isum,ii,jj

      if(inc.eq.1 .and. jnc.eq.1) then
         do j=1,jml
         do i=1,iml
            ii=lw(ip_x)+i+ism-5
            jj=ls(ip_y)+j+jsm-5
            pa1(i,j)=pa0(ii,jj)*texom(ii,jj,1)
     $           +(1.d0-texom(ii,jj,1))*999.d0

         enddo
         enddo

         return
      endif

      big=0.9d35
      sml=-0.9d35
      nx = im
      ny = jm

      isum=isumt(1)

      do n=1,isum
         xp(n) = slon0+(ipt(n,1)-1)*dxdeg0
         yp(n) = slat0+(jpt(n,1)-1)*dydeg0
         zp(n) = pa0(ipt(n,1),jpt(n,1))
      enddo

      do n = isum+1, ndpt
         xp(n) = xp(isum)
         yp(n) = yp(isum)
         zp(n) = zp(isum)
      enddo

      if(isum.ne.0) then

         x1=slon+dxdeg*dble(lw(ip_x)-1)
         y1=slat+dydeg*dble(ls(ip_y)-1)

         call zgrid(sp2,xp,yp,zp,nx,ny,x1,y1,dxdeg,dydeg,
     $        cay,nrng,ndpt,0)
         call smooth(sp2,nx,ny,nsm)

         do j=1,ny
         do i=1,nx
            if(sp2(i,j).gt.big.or.sp2(i,j).lt.sml) sp2(i,j)=999.d0
         enddo
         enddo


         do j=1,jml
         do i=1,iml
           pa1(i,j) = sp2(i,j)
         enddo
         enddo

      else

         do j = 1, jm
         do i = 1, im
            pa1(i,j) = 999.d0
         enddo
         enddo

      endif

      return
      end  subroutine intplh


!=================== inerpolation for um ============================

      subroutine intplum(pa1,pa0)

      use param
      use mod_mpi
      implicit none
#include "common.h"

      real(8),intent(in) ::  pa0(im0,jm0)
      real(8),intent(out) :: pa1(im,jm)
      real(8) :: sp2(im,jm),x1,y1,big,sml,smd(im0,jm0)
      integer :: i,j,n,isum,nx,ny,ii,jj

      if(inc.eq.1 .and. jnc.eq.1) then
         do j=1,jml
         do i=1,iml
            ii=lw(ip_x)+i+ism-5
            jj=ls(ip_y)+j+jsm-5
            pa1(i,j)=pa0(ii,jj)*exom(ii,jj,1)
     $           +(1.d0-exom(ii,jj,1))*999.d0

         enddo
         enddo

         return
      endif

      big=0.9d35
      sml=-0.9d35
      nx=im
      ny=jm

      isum=isumu(1)

      do j=1,jm0
      do i=1,im0
         if(exnom(i,j)>0) then
            smd(i,j)=pa0(i,j)
         else
            smd(i,j)=big
         endif
      enddo
      enddo

!      call bcast_dble(smd,im0*jm0)
!      if(npresm>0) then
!         call smooth(smd,im0,jm0,npresm)
!      endif

      do n=1,isum
         xp(n) = slonu0+(ipu(n,1)-1)*dxdeg0
         yp(n) = slatu0+(jpu(n,1)-1)*dydeg0
         zp(n) = smd(ipu(n,1),jpu(n,1))/csom(jpu(n,1))
      enddo

      do n = isum+1, ndpt
         xp(n) = xp(isum)
         yp(n) = yp(isum)
         zp(n) = zp(isum)
      enddo

      if(isum.ne.0) then
          x1=slonu+dxdeg*dble(lw(ip_x)-1)
          y1=slatu+dydeg*dble(ls(ip_y)-1)
         call zgrid(sp2,xp,yp,zp,nx,ny,x1,y1,dxdeg,dydeg,
     $        cay,nrng,ndpt,0)
         call smooth(sp2,nx,ny,nsm)

         do j=1,ny
         do i=1,nx
           if(sp2(i,j).gt.big.or.sp2(i,j).lt.sml) sp2(i,j)=999.d0
         enddo
         enddo

         do j = 1, jml
         do i = 1, iml
         if(i.le.im-1.and.j.le.jm-1.and.sp2(i+1,j+1).ne.999.d0)then
            pa1(i,j) = sp2(i,j)*cs(j)
         else
            pa1(i,j) = sp2(i,j)
         endif
         enddo
         enddo

      else

         do j=1,jm
         do i=1,im
            pa1(i,j)=999.d0
         enddo
         enddo

      endif


      return
      end subroutine intplum

!=================== inerpolation for vm ============================

      subroutine intplvm(pa1,pa0)

      use param
      use mod_mpi
      implicit none
#include "common.h"

      real(8),intent(in) ::  pa0(im0,jm0)
      real(8),intent(out) :: pa1(im,jm)
      real(8) ::  sp2(im,jm),x1,y1,big,sml,smd(im0,jm0)
      integer :: i,j,nx,ny,n,isum,ii,jj

      if(inc.eq.1 .and. jnc.eq.1) then
         do j=1,jml
         do i=1,iml
            ii=lw(ip_x)+i+ism-5
            jj=ls(ip_y)+j+jsm-5
            pa1(i,j)=pa0(ii,jj)*exom(ii,jj,1)
     $           +(1.d0-exom(ii,jj,1))*999.d0

         enddo
         enddo

         return
      endif

      big=0.9d35
      sml=-0.9d35
      nx=im
      ny=jm

      isum=isumu(1)

      do n=1,isum
         xp(n) = slonu0+(ipu(n,1)-1)*dxdeg0
         yp(n) = slatu0+(jpu(n,1)-1)*dydeg0
         zp(n) = pa0(ipu(n,1),jpu(n,1))
      enddo

      do j=1,jm0
      do i=1,im0
         if(exnom(i,j)>0) then
            smd(i,j)=pa0(i,j)
         else
            smd(i,j)=big
         endif
      enddo
      enddo

!      if(npresm>0) then
!         call smooth(smd,im0,jm0,npresm)
!      endif

      do n = isum+1, ndpt
         xp(n) = xp(isum)
         yp(n) = yp(isum)
         zp(n) = zp(isum)
      enddo

      if(isum.ne.0) then

        x1=slonu+dxdeg*dble(lw(ip_x)-1)
        y1=slatu+dydeg*dble(ls(ip_y)-1)

         call zgrid(sp2,xp,yp,zp,nx,ny,x1,y1,dxdeg,dydeg,
     $        cay,nrng,ndpt,0)
         call smooth(sp2,nx,ny,nsm)

         do j=1,ny
         do i=1,nx
            if(sp2(i,j).gt.big.or.sp2(i,j).lt.sml) sp2(i,j)=999.d0
         enddo
         enddo


         do j=1,jml
         do i=1,iml
           pa1(i,j) = sp2(i,j)
         enddo
         enddo

      else

         do j = 1, jm
         do i = 1, im
            pa1(i,j) = 999.d0
         enddo
         enddo

      endif


      return
      end subroutine intplvm

!=============== interpolation for u,v ==========================

      subroutine intplu(pa1,pa0)

      use param
      use mod_mpi

      implicit none
#include "common.h"

      real(8),intent(in) ::  pa0(im0,jm0,0:km+1)
      real(8),intent(out) :: pa1(im,jm,0:km+1)
      real(8) :: sp2(im,jm),x1,y1,big,sml
      integer :: i,j,k,ii,jj,nx,ny,isum,n

      if(inc.eq.1 .and. jnc.eq.1) then
         do k=1,km
         do j=1,jml
         do i=1,iml
            ii=lw(ip_x)+i+ism-5
            jj=ls(ip_y)+j+jsm-5
            pa1(i,j,k)=pa0(ii,jj,k)*exom(ii,jj,k)
     $           +(1.d0-exom(ii,jj,k))*999.d0

         enddo
         enddo
         enddo

         return
      endif

      big=0.9d35
      sml=-0.9d35

      do 10 k = 1, km

        nx=im
        ny=jm

         isum=isumu(k)

         do n=1,isum
            xp(n) = slonu0+(ipu(n,k)-1)*dxdeg0
            yp(n) = slatu0+(jpu(n,k)-1)*dydeg0
            zp(n) = pa0(ipu(n,k),jpu(n,k),k)
         enddo

         do n = isum+1, ndpt
            xp(n) = xp(isum)
            yp(n) = yp(isum)
            zp(n) = zp(isum)
         enddo

         if(isum.ne.0) then

            x1=slonu+dxdeg*dble(lw(ip_x)-1)
            y1=slatu+dydeg*dble(ls(ip_y)-1)

#ifdef DEBUG1
      if(ip .eq. 0) then
      write(*,*) 'level',k,' start'
      endif
#endif
#ifdef DEBUG5
      if(ip .eq. 0) then
      write(*,*) 'level',k,' start'
      endif
#endif
            call zgrid(sp2,xp,yp,zp,nx,ny,x1,y1,dxdeg,dydeg,
     $           cay,nrng,ndpt,0)
            call smooth(sp2,nx,ny,nsm)

            do j=1,ny
            do i=1,nx
               if(sp2(i,j).gt.big.or.sp2(i,j).lt.sml) sp2(i,j)=999.d0
            enddo
            enddo

            do j=1,jml
            do i=1,iml
                pa1(i,j,k)=sp2(i,j)
            enddo
            enddo
         else

            do j = 1, jml
            do i = 1, iml
               pa1(i,j,k) = 999.d0
            enddo
            enddo

         endif

 10   continue

      return
      end subroutine intplu




!=================== inerpolation for wsx and wsy============================

      subroutine intplws(pa1,pa0)
      use param
      use mod_mpi
      implicit none

#include "common.h"

      real(4),intent(in) ::  pa0(im0,jm0,nday)
      real(4),intent(out) :: pa1(im,jm,nday)
      real(8) ::  sp2(im,jm),x1,y1,big,sml
      integer :: i,j,n,nx,ny,m,isum,iflag

      big=0.9d35
      sml=-0.9d35
      nx = im
      ny = jm

      isum=isumu(1)

      do 10 m=1,nday

      do n=1,isum
         xp(n) = slonu0+(ipu(n,1)-1)*dxdeg0
         yp(n) = slatu0+(jpu(n,1)-1)*dydeg0
         zp(n) = pa0(ipu(n,1),jpu(n,1),m)
      enddo

      do n = isum+1, ndpt
         xp(n) = xp(isum)
         yp(n) = yp(isum)
         zp(n) = zp(isum)
      enddo

      if(isum.ne.0) then
          x1=slonu+dxdeg*dble(lw(ip_x)-1)
          y1=slatu+dydeg*dble(ls(ip_y)-1)

!      if(m.ge.18) then
!        iflag = 1
!      else
         iflag = 0
!      endif

         call zgrid(sp2,xp,yp,zp,nx,ny,x1,y1,dxdeg,dydeg,
     $        cay,nrng,ndpt,iflag)
         call smooth(sp2,nx,ny,nsm)

         do j=1,ny
         do i=1,nx
            if(sp2(i,j).gt.big.or.sp2(i,j).lt.sml) sp2(i,j)=999.d0
         enddo
         enddo

         do j=1,jm
         do i=1,im
           pa1(i,j,m) =sngl(sp2(i,j))
         enddo
         enddo

      else

         do j = 1, jm
         do i = 1, im
            pa1(i,j,m) = 999.
         enddo
         enddo

      endif

 10   continue

      return
      end

!=================== inerpolation for surface flux =========================


      subroutine intplsf(pa1,pa0)
      use param
      use mod_mpi
      implicit none
#include "common.h"

      real(4),intent(in) ::  pa0(im0,jm0,nday)
      real(4),intent(out) :: pa1(im,jm,nday)
      real(8) :: sp2(im,jm),x1,y1,big,sml
      integer :: i,j,n,m,nx,ny,isum

      big=0.9d35
      sml=-0.9d35
      nx = im
      ny = jm

      isum=isumt(1)

      do m=1,nday

      do n=1,isum
         xp(n) = slon0+(ipt(n,1)-1)*dxdeg0
         yp(n) = slat0+(jpt(n,1)-1)*dydeg0
         zp(n) = pa0(ipt(n,1),jpt(n,1),m)
      enddo

      do n = isum+1, ndpt
         xp(n) = xp(isum)
         yp(n) = yp(isum)
         zp(n) = zp(isum)
      enddo

      if(isum.ne.0) then
        x1=slon+dxdeg*dble(lw(ip_x)-1)
        y1=slat+dydeg*dble(ls(ip_y)-1)
         call zgrid(sp2,xp,yp,zp,nx,ny,x1,y1,dxdeg,dydeg,
     $        cay,nrng,ndpt,0)
         call smooth(sp2,nx,ny,nsm)

         do j=1,ny
         do i=1,nx
            if(sp2(i,j).gt.big.or.sp2(i,j).lt.sml) sp2(i,j)=999.d0
         enddo
         enddo


         do j=1,jm
         do i=1,im
           pa1(i,j,m) = sngl(sp2(i,j))
         enddo
         enddo

      else

         do j = 1, jm
         do i = 1, im
            pa1(i,j,m) = 999.
         enddo
         enddo
      endif
#ifdef DEBUG1
       if(ip .eq. imaster ) write(*,*) 'flux;',m,' finished'
#endif
      enddo

      return
      end


c*****************************************************************
c*          A two dimensional spline interpolation subroutine
c*          adopted to RIAMOM in april 2003 by Hideyuki Kawamura
c*                     at Japan Atomic Energy Research Institute
c*****************************************************************
c
c***********************************************************************
      subroutine zgrid(sp1,xp,yp,zp,nx,ny,x1,y1,
     &     dxs,dys,cayin,nrng,ndpt,iflag)
c***********************************************************************

!      implicit double precision (a-h,o-z)
      implicit none
      real*8 sp1(nx,ny),xp(ndpt),yp(ndpt),zp(ndpt),dxs,dys,
     &     cayin,zpij(ndpt),x1,y1
      real*8 zmax,zmin,dzmax,dzrms,eps,big,cay,zrange,zbase,hrange
      real*8 derzm,zsum,relax,z00,wgt,zim,zimm,zip,zipp,zjm,zjmm,zjp
      real*8 dzs,x,y,zpxy,zw,ze,zs,zn,a,b,c,d,zxy,delz,delzm,root,dzmaxf
      real*8 dzrms8,tpy,rootgs,relaxn,zijn,abz,dzrmsp,zjpp
      integer knxt(ndpt),imnew(ny),npg,iflag,kflag,iter,nnew,jms,ims
      integer nrng,nx,ny,i,j,k,n,itmax,npt,kk,ndpt,jmnew
!      real*8 zmax,zmin,dzmax,dzrms

      n=ndpt
      itmax=100
      eps=0.002d0
      big=0.9d35
      cay=cayin

c     initialize(!)
c***********************************************************************

      do j=1,ny
      do i=1,nx
         sp1(i,j)=0.d0
      enddo
      enddo

      do j=1,ny
         imnew(j)=0
      enddo
      jmnew=0

c     get zbase which will make all up values positive by 20*(zmax-zmin)
c***********************************************************************

      zmin=zp(1)
      zmax=zp(1)

      do k=2,n
        zmax=dmax1(zmax,zp(k))
        zmin=dmin1(zmin,zp(k))
      enddo

      call max_all1(zmax,1)
      call min_all1(zmin,1)

      zrange=zmax-zmin
      zbase=zrange*20.d0-zmin
      hrange=dmin1(dxs*(nx-1),dys*(ny-1))
      derzm=2.d0*zrange/hrange

c     set pointer array knxt
c***********************************************************************

      do k = n,1,-1
      knxt(k)=0
      i=(xp(k)-x1)/dxs+1.5d0
      j=(yp(k)-y1)/dys+1.5d0
      if(i .ge. 1 .and. i .le. nx .and. j .ge. 1 .and. j. le. ny) then
          knxt(k)=n+1
          if(sp1(i,j) .gt. 0.) then
            knxt(k) = nint(sp1(i,j))
          endif
        sp1(i,j)=k
      endif
      enddo

#ifdef DEBUG1
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
       if(ip.eq.6 .and. iflag .eq. 1) then
         i = 1
         j = 1
         write(*,*) 'knxt',sp1(i,j)
       endif
#endif

c     affix each data point zp to its nearby grid point. take avg zp if
c     more than one zp nearby the grid point. add zbase and complement.
c***********************************************************************
      do j=1,ny
      do i=1,nx
          sp1(i,j)=-1.d35
      enddo
      enddo

      do k=1,n
      if(knxt(k) .gt. 0) then
        npt=0
        zsum=0.d0
        i=(xp(k)-x1)/dxs+1.5d0
        j=(yp(k)-y1)/dys+1.5d0
        kk=k

 70     npt=npt+1
        zsum=zsum+zp(kk)
        knxt(kk)=-knxt(kk)
        kk=-knxt(kk)
        if(kk .le. n) goto 70

        sp1(i,j)=-zsum/dble(npt)-zbase
      endif
      enddo

#ifdef DEBUG1
       if(ip.eq.6 .and. iflag .eq. 1) then
        write(*,*) 'init',zbase,zrange,zmax,zmin
        do j = 1,5
          write(*,*) j,'  sp1;',(sp1(i,j),i=1,5)
        enddo
      endif
#endif

c     initially set each unset grid point to value of nearest known pt.
c***********************************************************************
      do iter=1,nrng
! exchange boundary for parallelization
      call exch_bnd(sp1,nx,ny)

      nnew=0
      do j = 1,ny
        imnew(j) = 0
      enddo
      jmnew = 0

      do i = 1,nx
      do j = 1,ny
      if(sp1(i,j) .lt. -big) then
        if(j.ne.1 .and. jmnew .le. 0) then
          zijn=dabs(sp1(i,j-1))
          if(zijn .lt. big) then
            imnew(j)=1
            jmnew=1
            sp1(i,j)=zijn
            nnew=nnew+1
          endif
        endif
        if(j.ne.ny) then
          zijn=dabs(sp1(i,j+1))
          if(zijn .lt. big) then
            imnew(j)=1
            jmnew=1
            sp1(i,j)=zijn
            nnew=nnew+1
          endif
         endif

        if(i .ne. 1 .and. imnew(j) .le. 0) then
          zijn=dabs(sp1(i-1,j))
          if(zijn .lt. big) then
            imnew(j)=1
            jmnew=1
            sp1(i,j)=zijn
            nnew=nnew+1
          endif
        endif
        if(i .ne. nx) then
          zijn=dabs(sp1(i+1,j))
          if(zijn .lt. big) then
            imnew(j)=1
            jmnew=1
            sp1(i,j)=zijn
            nnew=nnew+1
          endif
        endif
      else
        imnew(j)=0
        jmnew=0
      endif
      enddo
      enddo

      call isum_all1(nnew,1)

#ifdef DEBUG1
       if(ip.eq.6 .and. iflag .eq. 1) then
        write(*,*) 'iter;',nnew
        do j = 1,5
          write(*,*) j,'  sp1;',(sp1(i,j),i=1,5)
        enddo
        write(*,*)
      endif
#endif
#ifdef DEBUG1
       if(ip.eq.6 .and. iflag .eq. 1) then
        i = 1
        j = 1
        write(*,*) 'iter;',iter,nnew,sp1(i,j)
      endif
#endif
      if(nnew .le. 0) goto 200

      enddo

 200  continue

      do i=1,nx
      do j=1,ny
        abz=dabs(sp1(i,j))
        if(abz .ge. big) then
          sp1(i,j)=abz
        endif
      enddo
      enddo

#ifdef DEBUG1
       if(ip.eq.6 .and. iflag .eq. 1) then
        i = 1
        j = 1
        write(*,*) 'set',sp1(i,j)
      endif
#endif

c     improve the non-data points by applying point over-relaxation
c     using the laplace-spline equation (carres method is used)
c***********************************************************************
      dzrmsp=zrange
      relax=1.0
      do 4000 iter=1,itmax
      dzrms=0.0
      dzmax=0.0
      npg=0

      do j=1,ny
      do i=1,nx

      z00=sp1(i,j)
      if(z00 .lt. big .and. z00.ge.0) then
         wgt=0.d0
         zsum=0.d0

         ims=0
         if(i.gt.1) then
            zim=dabs(sp1(i-1,j))
            if(zim.lt.big) then
               ims=1
               wgt=wgt+1.0
               zsum=zsum+zim
               if(i.gt.2) then
                  zimm=dabs(sp1(i-2,j))
                  if(zimm.lt.big) then
                     wgt=wgt+cay
                     zsum=zsum-cay*(zimm-2.d0*zim)
                  endif
               endif
            endif
         endif
         if(i.lt.nx) then
            zip=dabs(sp1(i+1,j))
            if(zip.lt.big) then
               wgt=wgt+1.d0
               zsum=zsum+zip
               if(ims.gt.0) then
                  wgt=wgt+4.d0*cay
                  zsum=zsum+2.d0*cay*(zim+zip)
               endif
               if(i.lt.nx-1) then
                  zipp=dabs(sp1(i+2,j))
                  if(zipp.lt.big) then
                     wgt=wgt+cay
                     zsum=zsum-cay*(zipp-2.0*zip)
                  endif
               endif
            endif
         endif
         jms=0
         if(j.gt.1) then
            zjm=dabs(sp1(i,j-1))
            if(zjm.lt.big) then
               jms=1
               wgt=wgt+1.d0
               zsum=zsum+zjm
               if(j.gt.2) then
                  zjmm=dabs(sp1(i,j-2))
                  if(zjmm.lt.big) then
                     wgt=wgt+cay
                     zsum=zsum-cay*(zjmm-2.0*zjm)
                  endif
               endif
            endif
         endif

         if(j.lt.ny) then
            zjp=dabs(sp1(i,j+1))
            if(zjp.lt.big) then
               wgt=wgt+1.d0
               zsum=zsum+zjp
               if(jms.gt.0) then
                  wgt=wgt+4.d0*cay
                  zsum=zsum+2.d0*cay*(zjm+zjp)
               endif
               if(j.lt.ny-1) then
                  zjpp=dabs(sp1(i,j+2))
                  if(zjpp.lt.big) then
                     wgt=wgt+cay
                     zsum=zsum-cay*(zjpp-2.d0*zjp)
                  endif
               endif
            endif
         endif

         dzs=zsum/wgt-z00
         npg=npg+1
         dzrms=dzrms+dzs*dzs
         dzmax=dmax1(dabs(dzs),dzmax)
         sp1(i,j)=z00+dzs*relax
      endif
      enddo
      enddo

      call max_all1(dzmax,1)
      call sum_all1(dzrms,1)
      call isum_all1(npg,1)

! exchange boundary for parallelization
      call exch_bnd(sp1,nx,ny)


c     shift data points zp progressively back to their proper places as
c     the shape of surface z becomes evident.
c***********************************************************************

      if(mod(iter,10) .eq. 0) then
      do k=1,n

      knxt(k)=iabs(knxt(k))
      if(knxt(k) .gt. 0) then
        x=(xp(k)-x1)/dxs
        i=x+1.5d0
        x=x+1.d0-i
        y=(yp(k)-y1)/dys
        j=y+1.5d0
        y=y+1.d0-j
        zpxy=zp(k)+zbase
        z00=dabs(sp1(i,j))

      kflag = 0
      if(i .gt. 1) then
        zw=dabs(sp1(i-1,j))
      else
!        zw=1.0d35
      kflag = kflag + 1
      endif

      if(i .lt. nx) then
        ze=dabs(sp1(i+1,j))
      else
!        ze=1.0d35
      kflag = kflag + 2
      endif

!      if(ze .ge. big) then
!        if(zw .ge. big) then
!          ze=z00
!          zw=z00
!        else
!          ze=2.d0*z00-zw
!        endif
!      else
!        if(zw .ge. big) then
!          zw=2.d0*z00-ze
!        endif
!      endif

      if(kflag .eq. 3) then
        ze=z00
        zw=z00
      elseif(kflag .eq. 2) then
        ze=2.d0*z00-zw
      elseif(kflag .eq. 1) then
        zw=2.d0*z00-ze
      endif

      kflag = 0
!      zs=1.0d35
      if(j .gt. 1) then
        zs=dabs(sp1(i,j-1))
      else
        kflag = kflag + 1
      endif

!      zn=1.0d35
      if(j .lt. ny) then
        zn=dabs(sp1(i,j+1))
      else
        kflag = kflag + 2
      endif

      if(kflag .eq. 3) then
        zn=z00
        zs=z00
      elseif(kflag .eq. 2) then
        zn=2.d0*z00-zs
      elseif(kflag .eq. 1) then
        zs=2.d0*z00-zn
      endif
!      if(zn .ge. big) then
!        if(zs .ge. big) then
!          zn=z00
!          zs=z00
!        else
!          zn=2.0*z00-zs
!        endif
!      else
!        if(zs .ge. big) then
!          zs=2.0*z00-zn
!        endif
!      endif

      a=(ze-zw)*0.5d0
      b=(zn-zs)*0.5d0
      c=(ze+zw)*0.5d0-z00
      d=(zn+zs)*0.5d0-z00
      zxy=z00+a*x+b*y+c*x*x+d*y*y
      delz=z00-zxy
      delzm=derzm*(dabs(x)*dxs+dabs(y)*dys)*0.8d0

!      if(delz .gt. delzm) then
!        delz=delzm
!      endif
!      if(delz .lt. -delzm) then
!        delz=-delzm
!      endif
      delz = dmin1( dmax1(delz,-delzm), delzm)

      zpij(k)=zpxy+delz

      endif
      enddo


      do k=1,n
      if(knxt(k) .gt. 0) then
        npt=0
        zsum=0.d0
        i=(xp(k)-x1)/dxs+1.5d0
        j=(yp(k)-y1)/dys+1.5d0
        kk=k

 3420   npt=npt+1
        zsum=zsum+zpij(kk)
        knxt(kk)=-knxt(kk)
        kk=-knxt(kk)

        if(kk .le. n)  goto 3420

        sp1(i,j)=-zsum/dble(npt)
      endif
      enddo

! exchange boundary for parallelization
      call exch_bnd(sp1,nx,ny)

#ifdef DEBUG1
      if(iflag .eq. 1) then
        write(*,*) 'check1, ip,',ip,' iter=',iter
      endif
#endif

      endif


c     test for convergence
c***********************************************************************

      if(npg.eq.0) go to 4010
      dzrms=dsqrt(dzrms/dble(npg))
      if(dzrmsp.eq.0.0) go to 4010
      root=dzrms/dzrmsp
      dzrmsp=dzrms
      dzmaxf=dzmax/zrange

      if(mod(iter,10) .eq. 2) then
        dzrms8=dzrms
      endif

      if(mod(iter,10) .eq. 0) then
        root=dsqrt(dsqrt(dsqrt(dzrms/dzrms8)))
        if(root-0.9999d0 .lt. 0) then
#ifdef DEBUG1
      if(iflag .eq. 1) then
        write(*,*) 'check2, ip,',ip,dzmaxf,root,eps,dzmaxf/(1.-root),
     &   dzmax,zrange,dzrms8,dzrmsp,npg
      endif
#endif
        if(dzmaxf/(1.d0-root)-eps .le. 0 ) goto 4010

c     improve the relaxation factor.
c***********************************************************************

          if(mod(iter,20) .eq. 0) then
            if(relax-1.d0-root .lt. 0) then
              tpy=(root+relax-1.d0)/relax
              rootgs=tpy*tpy/root
              relaxn=2.d0/(1.d0+dsqrt(1.d0-rootgs))
              if(iter .ne. 60) then
                relaxn=relaxn-0.25d0*(2.d0-relaxn)
              endif
              relax=dmax1(relax,relaxn)
            endif
          endif
        endif
      endif
#ifdef DEBUG1
      if(iflag .eq. 1) then
        write(*,*) 'check3, ip,',ip,' iter=',iter
      endif
#endif
 4000 continue
 4010 continue

c     remove zbase from array z and return.
c***********************************************************************

      do i=1,nx
      do j=1,ny
      if(sp1(i,j) .lt. big) then
        sp1(i,j)=dabs(sp1(i,j))-zbase
      endif
      enddo
      enddo

#ifdef DEBUG1
       if(ip.eq.6 .and. iflag .eq. 1) then
        i = 1
        j = 1
        write(*,*) 'output ',sp1(i,j)
        write(*,*)
      endif
#endif
#ifdef DEBUG5
      if(iflag .eq.1) then
         spmax=0.
         spmin=0.
         ispmax=0
         ispmin=0
         jspmax=0
         jspmin=0
         do i=2,nx-1
         do j=2,ny-1
            if(sp1(i,j) .lt. big) then
               if(sp1(i,j).gt.spmax) then
                  spmax=sp1(i,j)
                  ispmax=i-1
                  jspmax=j-1
               endif
               if(sp1(i,j).lt.spmin) then
                  spmin=sp1(i,j)
                  ispmin=i-1
                  jspmin=j-1
               endif
            endif
         enddo
         enddo
         write(*,*)
         write(*,*) 'output'
         write(*,'("max=",e8.2," (",i3,",",i3,")")' )
     $        spmax,ispmax,jspmax
         write(*,'("min=",e8.2," (",i3,",",i3,")")' )
     $        spmin,ispmin,jspmin
      endif
#endif
      return
      end

c***********************************************************************
      subroutine smooth(sp1,nx,ny,nsm)
c***********************************************************************

      implicit none
!      implicit double precision(a-h,o-z)
      real*8 sp1(nx,ny),big,r,zij,del2,del2x,del2y
      integer i,j,it
      integer,intent(in) :: nx,ny,nsm

      if(nsm .gt. 0) then
        big=0.9d35
        r=0.25d0/4.0d0

!        call exch_bnd(sp1,nx,ny)

        do it=1,nsm

          do j=2,ny-1
          do i=2,nx-1
            zij=sp1(i,j)

            if(zij .lt. big) then
              del2=0.0
                del2x=sp1(i-1,j)+sp1(i+1,j)-zij-zij
                if(del2x .lt. big) then
                  del2=del2x
                endif

                del2y=sp1(i,j-1)+sp1(i,j+1)-zij-zij
                if(del2y .lt. big) then
                  del2=del2+del2y
                endif

              sp1(i,j)=zij+del2*r

            endif
          enddo
          enddo

          call exch_bnd(sp1,nx,ny)

        enddo
      endif

#ifdef DEBUG1
       if(ip.eq.6) then
        i = 1
        j = 1
        write(*,*) 'smooth ',sp1(i,j)
        write(*,*)
      endif
#endif


      return
      end



!================== count side bounday data number ======================
      subroutine sbnum
      use param
      use mod_time
      use mod_mpi
      implicit none
#include "common.h"
      real(8) :: chour

      if(obc_clim==1) then
         chour=dble(nsec)/36.d2+idint(ahourbd(1)/nday/24)*nday*24.d0
         if(chour<ahourbd(1)) then
            chour=chour+dble(nday*24)
            nrct=1
         endif
         if(chour<ahourbd(nrct) .or. nrct+1>nrbd) nrct=1
         do while(chour>ahourbd(nrct+1))
            nrct=nrct+1

            if (nrct+1>nrbd)then
               chour=chour-dble(nday*24)
               nrct=1
            endif
         enddo
      else
         chour=ahour
         do while(chour>ahourbd(nrct+1))
            nrct=nrct+1
            if (nrct+1>nrbd)then
               write(*,*) 'STOP : Side boundary data are short.'
               call wrcontin
               stop
            endif
         enddo
      endif
#ifndef ASSIM
      if(ip==imaster .and. mod(nkai,nxday).eq.0) then
         write(*,*) 'chour=',chour,
     &        ' ahourbd(nrct)=',ahourbd(nrct),
     &        ' ahourbd(nrct+1)=',ahourbd(nrct+1)
         write(*,*) 'nrct=',nrct,'nkai=',nkai,
     &        'nsec=',nsec
      endif
#endif

!#ifdef ASSIM
!      if(chour<ahourbd(nrct)*3.6d3) then
!         nrct=nrct-1
!         goto 51
!      endif
!#endif

      mb1=nrct
      mb2=mb1+1

      cbf=(chour-ahourbd(mb1))/(ahourbd(mb2)-ahourbd(mb1))
#ifdef DEBUG
      if(cbf<0 .or. cbf>1) then
         if(ip==imaster) write(*,*) 'err : cbf=',cbf
         stop
      endif
#endif

      return
      end

!--------------------------------------------------------------------
      subroutine read_fd_nest(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,n
      real arr(im,jm,nday)
      real buf_in(img,jmg,nday)
c
      if(ip.eq.imaster) then
        read(nunit) buf_in
        write(*,*) 'read buf_in'
      endif
c
      call bcast_real(buf_in,img*jmg*nday)
c
      do n = 1,nday
      do j = 1,jml
      do i = 1,iml
         arr(i,j,n) = buf_in(lw(ip_x)-1+i,ls(ip_y)-1+j,n)
      enddo
      enddo
      enddo

      return
      end
!------------------------------------------------------------------------
      subroutine write_fd_nest(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,k,m,ii,jj,ipxx,ipyy

      real arr(im,jm,nday)
      real sbuf((im-4)*(jm-4)*nday),rbuf((im-4)*(jm-4)*nday*np)
      real x3d8(img,jmg,nday)

      if(ip.eq.imaster) then
        do k = 1,nday
        do j = 1,jmg
        do i = 1,img
          x3d8(i,j,k) = 0.
        enddo
        enddo
        enddo
      endif
c
#ifdef NO_MPI
      write(nunit) arr
#else
c
      do k = 1,nday
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*nday*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*nday,mpi_real,
     &     rbuf,(im-4)*(jm-4)*nday,mpi_real,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
        do k = 1,nday
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(k-1)+(im-4)*nday*(j-3)
     &        +(im-4)*(jm-4)*nday*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
           if(ii.le.img .and. jj.le.jmg) then
             x3d8(ii,jj,k) = rbuf(n)
           endif
        enddo
        enddo
        enddo
      enddo
        write(nunit) x3d8
!        write(*,*) 'write x3d8'
      endif
#endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      return
      end

!------------------------------------------------------
      subroutine exch_bnd(sp1,nx,ny)
      use param
      use mod_mpi
      implicit none

      integer nx,ny,i,j
      real*8 sp1(nx,ny)
      real*8 bufx(2,ny),bufy(nx,2)

      if(nx .ne. im .and. ny. ne. jm) then
        write(*,*) ' cannot exchange boundary values in nest.F'
        stop
      endif

      if(ip_x .eq. 0) then
        do j = 1,ny
        do i = 1,2
          bufx(i,j) = sp1(i,j)
        enddo
        enddo
      endif
      if(ip_x .eq. ipe-1) then
        do j = 1,ny
        do i = 1,2
          bufx(i,j) = sp1(iml-i+1,j)
        enddo
        enddo
      endif
      if(ip_y .eq. 0) then
        do j = 1,2
        do i = 1,nx
          bufy(i,j) = sp1(i,j)
        enddo
        enddo
      endif
      if(ip_y .eq. jpe-1) then
        do j = 1,2
        do i = 1,nx
          bufy(i,j) = sp1(i,jml+1-j)
        enddo
        enddo
      endif

      call exch_t2d_s1(sp1,1)
      call exch_t2d_s2(sp1,2)
      call exch_t2d_n1(sp1,3)
      call exch_t2d_n2(sp1,4)

      call wait_t2d_s1(sp1,1)
      call wait_t2d_s2(sp1,2)
      call wait_t2d_n1(sp1,3)
      call wait_t2d_n2(sp1,4)

      call exch_t2d_e1(sp1,5)
      call exch_t2d_e2(sp1,6)
      call exch_t2d_w1(sp1,7)
      call exch_t2d_w2(sp1,8)

      call wait_t2d_e1(sp1,5)
      call wait_t2d_e2(sp1,6)
      call wait_t2d_w1(sp1,7)
      call wait_t2d_w2(sp1,8)

      if(ip_x .eq. 0) then
        do j = 1,ny
        do i = 1,2
          sp1(i,j) = bufx(i,j)
        enddo
        enddo
      endif
      if(ip_x .eq. ipe-1) then
        do j = 1,ny
        do i = 1,2
          sp1(iml-i+1,j) =  bufx(i,j)
        enddo
        enddo
      endif
      if(ip_y .eq. 0) then
        do j = 1,2
        do i = 1,nx
          sp1(i,j) = bufy(i,j)
        enddo
        enddo
      endif
      if(ip_y .eq. jpe-1) then
        do j = 1,2
        do i = 1,nx
          sp1(i,jml+1-j) = bufy(i,j)
        enddo
        enddo
      endif


      return
      end


!---- end nest.F
#endif