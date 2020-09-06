ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     subroutine surf_bulk                                          c
c                                                                   c
c     $Id: surf_bulk.F90 11 2008-12-11 05:08:04Z ishikawa $
c                                                                   c
c     making surface flux from bulk formulae                        c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine surf_bulk

      use param
      use mod_mpi
      implicit none
#include "common.h"

#ifndef EQ100M

      real*8 fca,fcb,stdev_wind(im,jm),dpt_t2m(im,jm),dconduc(im,jm),
     &     temp_2m(im,jm),senh_i,
     &     tt,e_w_d(im,jm),q_a(im,jm),e_w,q_w,e_i,
     &     q_i,lath_i,lwr_o(im,jm),lwr_i,swr_i(im,jm),
     &     conduc(im,jm),win(im,jm),h(im,jm),dlwr_i(im,jm),factw,
     &     dsenh_i,dnetq_i,de_i,dq_i,dlath_i,t_old(im,jm),
     &     delta_ttop,
#ifdef NCEPSF

#ifdef DLWRF
     &     dswrf(im,jm),dlwrf(im,jm)
#else
     &     dswrf(im,jm),cloud(im,jm)
#endif

#else
     &     cloud(im,jm),tot_sol(im,jm)
#endif

!      integer*8 nkai_year,nsec,ndayb,ndaya
#ifdef ROKKA
      integer,save :: nhra=2,nhrb=1
#if defined(EXP2006) || defined(EXP2007)
      integer,save :: nwrsa=2,nwrsb=1
      real(8) :: fca_ws,fcb_ws
#endif
#else
      integer*8 nkai_year,ndayb,ndaya
#endif
      integer it,maxit,i,j
#ifdef EVAPO
      real*8 evapo
#ifndef NESTED
      real*8 dum(2),totarea,htmnpice,tscalewf
#endif
#endif


!for assim
      real*8 alp_m,alp_w,alp_h
      common /cntl1/ alp_m(im,jm),alp_w(im,jm),alp_h(im,jm)
      real*8 lath_o,senh_o
      common /surf/ lath_o(im,jm),senh_o(im,jm)


#ifndef ASSIM
      do j = 1,jm
      do i = 1,im
        alp_m(i,j) = 1.0
        alp_w(i,j) = 1.0
        alp_h(i,j) = 1.0
      enddo
      enddo
#endif

!   initialize flux
      do j = 1,jm
      do i = 1,im
         wflux(i,j)=0.d0
         swr_o(i,j)=0.d0
         netq_o(i,j)=0.d0
         netq_i(i,j)=0.d0
         fricv3(i,j)=0.d0
         wsx(i,j) = 0.d0
         wsy(i,j) =0.d0
      enddo
      enddo
#ifndef ROKKA

!   this is daily forcing case

#ifndef NESTED
      nkai_year=mod(nkai-1,nxyear)+1
      nsec_surf=mod((nkai_year-1)*nint(dtuv),nday*24*60*60)
#endif
#ifdef DEBUG11
      if(mod(nkai,nxday).eq.1 .and. matsn2.ne.1) then
          if(ip.eq.imaster) then
            write(*,'("Julian day = ",i3)') nsec/86400+1
         endif
      endif
#endif
      ndayb=(nsec+12*60*60)/(24*60*60)
      ndaya=ndayb+1
      fca=dble(nsec-ndayb*24*60*60+12*60*60)/(24.*60.*60.)
      fcb=1.-fca
#ifdef CLIMAT
      if(ndaya==nday+1)ndaya=1
      if(ndayb==0)ndayb=nday
#endif
#else
#if !defined(EXP2006) && !defined(EXP2007)
!   6-hourly (0:00, 6:00, 12:00, 18:00)
         nhrb=nsec/60/60/6+1
         nhra=nhrb+1
         fca=dble(nsec-(nhrb-1)*6*60*60)/(6.d0*60.d0*60.d0)
         fcb=1.d0-fca

!   3-hourly (0:00 3:00 6:00 9:00 12:00 15:00 18:00 21:00) for wind stress
!    MSM is updatd after Mar 2006
#else
!  for scolar data
#ifdef DEBUG13
         if(ip==imaster)
     $   write(*,
     $     "('nsec_scf(nhrb)=',i,2x,'nsec=',i,2x,'nsec_scf(nhra)=',i)")
     $        nsec_scf(nhrb),nsec,nsec_scf(nhra)
#endif
         if(nsec_scf(nhrb)>nsec) then
            nhra=2
            nhrb=1
         endif

         do while(nsec_scf(nhra)<nsec)
            nhrb=nhrb+1
            nhra=nhrb+1
#ifdef DEBUG13
            if(ip==imaster) then
               print *,'nhra=',nhra
               print *,'nsec=',nsec,'nsec_scf=',nsec_scf(nhra)
            endif
#endif
#ifdef DEBUG
            if(nhra>nsf*3) then
               if(ip==imaster) print *,"surface flux data are short"
               stop
            endif
#endif
         enddo


         fca=dble(nsec-nsec_scf(nhrb))
     $        /dble(nsec_scf(nhra)-nsec_scf(nhrb))
         fcb=1.d0-fca

! for wind stress data

#ifdef DEBUG13
         if(ip==imaster)
     $   write(*,
     $    "('nsec_scf(nwrsb)=',i,2x,'nsec=',i,2x,'nsec_scf(nwrsa)=',i)")
     $      nsec_scf(nwrsb),nsec,nsec_scf(nwrsa)
#endif

         if(nsec_ws(nwrsb)>nsec) then
            nwrsa=2
            nwrsb=1
         endif

         do while(nsec_ws(nwrsa)<nsec)
            nwrsb=nwrsb+1
            nwrsa=nwrsb+1
#ifdef DEBUG13
            if(ip==imaster) then
               print *,'nwrsa=',nwrsa
               print *,'nsec=',nsec,'nsec_ws=',nsec_ws(nwrsa)
            endif

#endif
#ifdef DEBUG
            if(nwrsa>nws*3) then
               if(ip==imaster) print *,"wind stress data are short"
               stop
            endif
#endif
         enddo

         fca_ws=dble(nsec-nsec_ws(nwrsb))
     $        /dble(nsec_ws(nwrsa)-nsec_ws(nwrsb))
         fcb_ws=1.d0-fca_ws
#ifdef DEBUG13
         if(ip==imaster) then
            write(*,*) 'nsec=',nsec
            write(*,*) 'nhra=',nhra,' nhrb=',nhrb,
     $           ' fca=',fca,' fcb=',fcb
            write(*,*) 'nwrsa=',nwrsa
            write(*,*) 'nsec_ws(nwrsa)=',nsec_ws(nwrsa)
            write(*,*) 'fca_ws=',fca_ws,' fcb_ws=',fcb_ws
         endif
#endif
#endif
#ifdef CLIMAT
         if(nhra>nday*4) nhra=1
#endif
#endif

      do j=3,jml-2
      do i=3,iml-2

!   interpolate atmospheric condition
#ifndef ROKKA
         win(i,j)= tex(i,j,1)*dabs
     &     (fcb*sc_wind_d(i,j,ndayb)+fca*sc_wind_d(i,j,ndaya))/100.
         stdev_wind(i,j)=tex(i,j,1)*dabs
     &        (fcb*stdev_w_d(i,j,ndayb)+fca*stdev_w_d(i,j,ndaya))
         dpt_t2m(i,j)=tex(i,j,1)*
     &        (fcb*dpt_t2m_d(i,j,ndayb)+fca*dpt_t2m_d(i,j,ndaya))
         temp_2m(i,j)=tex(i,j,1)*
     &        (fcb*temp_2m_d(i,j,ndayb)+fca*temp_2m_d(i,j,ndaya))
#ifdef NCEPSF
         dswrf(i,j)=tex(i,j,1)*
     &        (fcb*dswrf_d(i,j,ndayb)+fca*dswrf_d(i,j,ndaya))
#ifdef DLWRF
         dlwrf(i,j)=tex(i,j,1)*
     &        (fcb*dlwrf_d(i,j,ndayb)+fca*dlwrf_d(i,j,ndaya))
#else
         cloud(i,j)=tex(i,j,1)*
     &        (fcb*cloud_d(i,j,ndayb)+fca*cloud_d(i,j,ndaya))
#endif
#else
         cloud(i,j)=tex(i,j,1)*
     &        (fcb*cloud_d(i,j,ndayb)+fca*cloud_d(i,j,ndaya))
         tot_sol(i,j)=tex(i,j,1)*
     &        (fcb*tot_sol_d(i,j,ndayb)+fca*tot_sol_d(i,j,ndaya))
#endif
#ifdef DEBUG1
      if(nkai.ge.1670 .and. nkai .le. 1677 .and. ip.eq.14
     & .and. i.eq.16 .and. j.eq.25) then
        write(*,*) fcb,fca,ndayb,ndaya,ntday_year,nday,nsec
        write(*,*) nkai,'win',win(i,j),sc_wind_d(i,j,ndayb),
     &    sc_wind_d(i,j,ndaya)
        write(*,*) nkai,'dev',stdev_wind(i,j),stdev_w_d(i,j,ndayb),
     &    stdev_w_d(i,j,ndaya)
      endif
#endif
#ifdef DEBUGS
         if(nkai.ge.nn-10.and.nkai.le.nn+10.and.
     &       lw(ip_x)+i-1.eq.ii.and.ls(ip_y)+j-1.eq.jj)then
            write(*,*)' '
            write(*,*)'omip(wind/sd_win/dpt_t2m/t2m/cloud/tot_sol)'
            write(*,*)win(i,j),stdev_wind(i,j),dpt_t2m(i,j),
     &           temp_2m(i,j),cloud(i,j),tot_sol(i,j)
         endif
#endif

#else
! for ROKKA
         win(i,j)= tex(i,j,1)*dabs
     &        (fcb*sc_wind_d(i,j,nhrb)+fca*sc_wind_d(i,j,nhra))/100.
         stdev_wind(i,j)=tex(i,j,1)*dabs
     &        (fcb*stdev_w_d(i,j,nhrb)+fca*stdev_w_d(i,j,nhra))
         dpt_t2m(i,j)=tex(i,j,1)*
     &        (fcb*dpt_t2m_d(i,j,nhrb)+fca*dpt_t2m_d(i,j,nhra))
         temp_2m(i,j)=tex(i,j,1)*
     &           (fcb*temp_2m_d(i,j,nhrb)+fca*temp_2m_d(i,j,nhra))
#ifdef NCEPSF
         dswrf(i,j)=tex(i,j,1)*
     &        (fcb*dswrf_d(i,j,nhrb)+fca*dswrf_d(i,j,nhra))
#ifdef DLWRF
         dlwrf(i,j)=tex(i,j,1)*
     &        (fcb*dlwrf_d(i,j,nhrb)+fca*dlwrf_d(i,j,nhra))
#else
         cloud(i,j)=tex(i,j,1)*
     &        (fcb*cloud_d(i,j,nhrb)+fca*cloud_d(i,j,nhra))
#endif
#else
         cloud(i,j)=tex(i,j,1)*
     &        (fcb*cloud_d(i,j,nhrb)+fca*cloud_d(i,j,nhra))
         tot_sol(i,j)=tex(i,j,1)*
     &        (fcb*tot_sol_d(i,j,nhrb)+fca*tot_sol_d(i,j,nhra))
#endif
#endif

!   wind stress
!   in case under ice; modified in sub. "ice-pre"
!   spinup gradually

#ifdef ROKKA
#if defined(EXP2006) || defined(EXP2007)
         wsx(i,j) = alp_m(i,j)*factw*ex(i,j,1)
     &        *(fcb_ws*wsx_d(i,j,nwrsb)+fca_ws*wsx_d(i,j,nwrsa))
         wsy(i,j) = alp_m(i,j)*factw*ex(i,j,1)
     &        *(fcb_ws*wsy_d(i,j,nwrsb)+fca_ws*wsy_d(i,j,nwrsa))
#else
         wsx(i,j) = alp_m(i,j)*factw*ex(i,j,1)
     &        *(fcb*wsx_d(i,j,nhrb)+fca*wsx_d(i,j,nhra))
         wsy(i,j) = alp_m(i,j)*factw*ex(i,j,1)
     &        *(fcb*wsy_d(i,j,nhrb)+fca*wsy_d(i,j,nhra))
#endif
         wflux(i,j)=fcb*wflux_d(i,j,nhrb)+fca*wflux_d(i,j,nhra)
#else
         wsx(i,j) = alp_m(i,j)*factw*ex(i,j,1)
     &        *(fcb*wsx_d(i,j,ndayb)+fca*wsx_d(i,j,ndaya))
         wsy(i,j) = alp_m(i,j)*factw*ex(i,j,1)
     &        *(fcb*wsy_d(i,j,ndayb)+fca*wsy_d(i,j,ndaya))

!   fresh water flux is not bulk type now
         wflux(i,j)=fcb*wflux_d(i,j,ndayb)+fca*wflux_d(i,j,ndaya)
#endif

      enddo
      enddo

#ifdef ICE
!   determine the ice top temperature

      do j=3,jml-2
      do i=3,iml-2
      if(tex(i,j,1)*aice(i,j).gt.0.d0
     &        .and.temp_2m(i,j)-273.15.lt.t(i,j,1))then

        e_w_d(i,j)=a_w*exp((b_w-(dpt_t2m(i,j)-273.15)/d_w)
     &    *(dpt_t2m(i,j)-273.15)/((dpt_t2m(i,j)-273.15)+c_w))
        q_a(i,j)=0.622*e_w_d(i,j)/(p_a-0.378*e_w_d(i,j))

#ifdef NCEPSF
        swr_i(i,j)=aice(i,j)*dswrf(i,j)*(1.-albedo_i)*1.e4
#else
        swr_i(i,j)=aice(i,j)*tot_sol(i,j)*1.e4
     &       *(1.-0.62*cloud(i,j))*(1.-albedo_i)
#endif

        dlwr_i(i,j)=-aice(i,j)*4.*eps_w*cboltz*temp_2m(i,j)**3

        h(i,j)=volice(i,j)/(areat(i,j,1)*aice(i,j))*1.d-2
        if(h(i,j).gt.4.0d0) then
          h(i,j)=4.d0+(h(i,j)-4.0d0)*26.d0
        endif
        dconduc(i,j)=-aice(i,j)*k_i/h(i,j)
c
ccc   iteration for t_top
        t_old(i,j) = t_top(i,j)
        if(t_top(i,j).eq.t(i,j,1))
     &       t_top(i,j)=(temp_2m(i,j)-273.15+t(i,j,1))/2.
      else
        t_top(i,j)=t(i,j,1)*tex(i,j,1)
        h(i,j)=0.
        t_old(i,j)=t_top(i,j)
      endif
c
#ifdef DEBUGI
      if(nkai.ge.nn-10.and.nkai.le.nn+10.and.
     &     i.eq.ii.and.ls(ip)+j-1.eq.jj)then
         write(*,*)' '
         write(*,*)nkai,i,ls(ip)+j-1
      endif
#endif
      enddo
      enddo

      do j=3,jml-2
      do i=3,iml-2

      maxit=0

 998  continue
      it=0
      if(tex(i,j,1)*aice(i,j).gt.0.d0
     &        .and.temp_2m(i,j)-273.15.lt.t(i,j,1))then

        e_i=a_i*exp((b_i-t_top(i,j)/d_i)*t_top(i,j)/(t_top(i,j)+c_i))
        q_i=0.9815*0.622*e_i/(p_a-0.378*e_i)
        lath_i=aice(i,j)*rho_a*l_w*c_la_w*win(i,j)*(q_a(i,j)-q_i)
        de_i=e_i*(b_i-t_top(i,j)/d_i)*t_top(i,j)/(t_top(i,j)+c_i)
     &       *( (b_i-t_top(i,j)/d_i*2.d0)/(t_top(i,j)+c_i)
     &       -(b_i-t_top(i,j)/d_i)*t_top(i,j)/(t_top(i,j)+c_i)**2 )
        dq_i=0.9815*0.622*de_i*p_a/(p_a-0.378*de_i)**2
        dlath_i=aice(i,j)*rho_a*l_w*c_la_w*win(i,j)*(-dq_i)
        if(lath_i.gt.0.d0)then
           lath_i=0.d0
           dlath_i=0.d0
        endif
c
#ifdef DLWRF
        lwr_i=aice(i,j)*( dlwrf(i,j)
     &       -4.*eps_w*cboltz*temp_2m(i,j)**3
     &       *(t_top(i,j)+273.15-temp_2m(i,j)) )
#else
        lwr_i=-aice(i,j)*( eps_w*cboltz*temp_2m(i,j)**4
     &       *(0.39-0.05*sqrt(e_w_d(i,j)))*(1.-xkai(j)*cloud(i,j)**2)
     &       +4.*eps_w*cboltz*temp_2m(i,j)**3
     &       *(t_top(i,j)+273.15-temp_2m(i,j)) )
#endif
        conduc(i,j)=aice(i,j)*k_i*(t(i,j,1)-t_top(i,j))/h(i,j)
c
        tt=temp_2m(i,j)-(t_top(i,j)+273.15)
        if(t_top(i,j).lt.temp_2m(i,j)-273.15d0)then
           senh_i=aice(i,j)*rho_a*c_p_a*c_stab*win(i,j)
     &        *dmax1(tt-dtskin,0.d0)
           dsenh_i=-aice(i,j)*rho_a*c_p_a*c_stab*win(i,j)
        else
           senh_i=aice(i,j)*rho_a*c_p_a*c_unst*win(i,j)
     &      *dmin1(tt+dtskin,0.d0)
           dsenh_i=-aice(i,j)*rho_a*c_p_a*c_unst*win(i,j)
        endif
c
        netq_i(i,j)=swr_i(i,j)+lwr_i+senh_i+lath_i+conduc(i,j)
        if(netq_i(i,j).gt.0.d0) then
#ifdef NCEPSF
          swr_i(i,j)=aice(i,j)*dswrf(i,j)*(1.-albedo_melt)*1.d4
#else
          swr_i(i,j)=aice(i,j)*tot_sol(i,j)*1.d4
     &       *(1.-0.62*cloud(i,j))*(1.-albedo_melt)
#endif
        endif
c
        netq_i(i,j)=swr_i(i,j)+lwr_i+senh_i+lath_i+conduc(i,j)
        dnetq_i=dlwr_i(i,j)+dsenh_i+dlath_i+dconduc(i,j)
        delta_ttop=-netq_i(i,j)/dnetq_i
        t_top(i,j)=t_top(i,j)+delta_ttop
        t_top(i,j)=dmax1(t_top(i,j),temp_2m(i,j)-273.15d0)
        t_top(i,j)=dmin1(t_top(i,j),t(i,j,1))
c
      if(dabs(t_top(i,j)-t_old(i,j)).ge.t_ice_lim
     &       .or.maxit.le.10)then
        t_old(i,j) = t_top(i,j)
        it=1
      endif
c
      endif
c
      if(it.ne.0)then
         maxit=maxit+1
         if(maxit.le.100)then
            goto 998
         endif
      endif
c
#ifdef DEBUGI
         if(tex(i,j,1)*aice(i,j).gt.0.d0.and.
     &        temp_2m(i,j)-273.15.lt.t(i,j,1).and.
     &        temp_2m(i,j)-273.15.lt.-10..and.
     &        (temp_2m(i,j)-273.15.eq.t_top(i,j)
     &        .or.t(i,j,1).eq.t_top(i,j)))then
            write(*,*)' ta/tt/to: ',
     &           nkai,lw(ip_x)+i-1,ls(ip_y)+j-1,
     &           temp_2m(i,j)-273.15,t_top(i,j),t(i,j,1),
     &           netq_i(i,j),dnetq_i,delta_ttop,maxit,it
            write(*,*)'                        ',
     &           swr_i(i,j),lwr_i,senh_i,lath_i,conduc(i,j)
            write(*,*)'                        ',
     &           0.,dlwr_i(i,j),dsenh_i,dlath_i,dconduc(i,j)
            write(*,*)'            senh..            ',
     &           aice(i,j),rho_a,c_p_a,c_unst,c_cnst,win(i,j),dtskin
         endif
#endif
      enddo
      enddo
#endif
c
c
ccc   calculate surface forcing
c
      do j=3,jml-2
      do i=3,iml-2
      if(tex(i,j,1).eq.1.d0) then
ccc   sensible heat for ocean
        tt=temp_2m(i,j)-(t(i,j,1)+273.15)
        if(tt.ge.0.d0)then
           tt=dmax1(tt-dtskin,0.d0)
        else
           tt=dmin1(tt+dtskin,0.d0)
        endif
        senh_o(i,j)=1.d-4*rho_a*c_p_a*c_stab*win(i,j)*tt
c
#ifdef ICE
ccc   sensible heat for ice
        tt=temp_2m(i,j)-(t_top(i,j)+273.15)
        if(t_top(i,j).lt.temp_2m(i,j)-273.15d0)then
           tt=dmax1(tt-dtskin,0.d0)
           senh_i=aice(i,j)*1.d-4*rho_a*c_p_a*c_stab*win(i,j)*tt
        else
           tt=dmin1(tt+dtskin,0.d0)
           senh_i=aice(i,j)*1.d-4*rho_a*c_p_a*c_unst*win(i,j)*tt
        endif
#endif
c
ccc   latent heat from atmos. vapor presure
        e_w_d(i,j)=a_w*exp((b_w-(dpt_t2m(i,j)-273.15)/d_w)
     &  *(dpt_t2m(i,j)-273.15)/((dpt_t2m(i,j)-273.15)+c_w))
        q_a(i,j)=0.622*e_w_d(i,j)/(p_a-0.378*e_w_d(i,j))
c
ccc   for ocean
        e_w=a_w*exp((b_w-t(i,j,1)/d_w)
     &   *t(i,j,1)/(t(i,j,1)+c_w))
        q_w=0.9815*0.622*e_w/(p_a-0.378*e_w)
        lath_o(i,j)=1.d-4*rho_a*l_w*c_la_w*win(i,j)*(q_a(i,j)-q_w)
        lath_o(i,j)=dmin1(lath_o(i,j),0.d0)
c
#ifdef ICE
ccc   for ice
        e_i=a_i*exp((b_i-t_top(i,j)/d_i)
     &   *t_top(i,j)/(t_top(i,j)+c_i))
        q_i=0.9815*0.622*e_i/(p_a-0.378*e_i)
        lath_i=aice(i,j)*1.d-4*rho_a*l_w*c_la_w*win(i,j)
     &     *(q_a(i,j)-q_i)
        lath_i=dmin1(lath_i,0.d0)
#endif

!   longwave radiation for ocean
#ifdef DLWRF
        lwr_o(i,j)=dlwrf(i,j)
     &       -1.d-4*4.*eps_w*cboltz*temp_2m(i,j)**3*(t(i,j,1)+273.15-
     &       temp_2m(i,j))
#else
        lwr_o(i,j)=-1.d-4*( eps_w*cboltz*temp_2m(i,j)**4
     &       *(0.39-0.05*sqrt(e_w_d(i,j)))
     &          *(1.-xkai(j)*cloud(i,j)**2)
     &       +4.*eps_w*cboltz*temp_2m(i,j)**3
     &       *(t(i,j,1)+273.15-temp_2m(i,j)) )
#endif

#ifdef ICE
!   longwave radiation for ice
#ifdef DLWRF
        lwr_i=aice(i,j)*(dlwrf(i,j)
     &       -1.d-4*4.*eps_w*cboltz*temp_2m(i,j)**3
     &       *(t_top(i,j)+273.15-temp_2m(i,j)) )
#else
        lwr_i=-aice(i,j)*1.d-4*( eps_w*cboltz*temp_2m(i,j)**4
     &       *(0.39-0.05*sqrt(e_w_d(i,j)))*(1.-xkai(j)*cloud(i,j)**2)
     &       +4.*eps_w*cboltz*temp_2m(i,j)**3
     &       *(t_top(i,j)+273.15-temp_2m(i,j)) )
#endif
#endif

!   short wave radiation for ocean
#ifdef ICE
#ifdef NCEPSF
        swr_o(i,j)=(1.-aice(i,j))*dswrf(i,j)*(1.-albedo_w)
        swr_i(i,j)=aice(i,j)*dswrf(i,j)*(1.-albedo_i)
#else
        swr_o(i,j)=(1.-aice(i,j))*tot_sol(i,j)
     &       *(1.-0.62*cloud(i,j))*(1.-albedo_w)
!   for ice
        swr_i(i,j)=aice(i,j)*tot_sol(i,j)
     &       *(1.-0.62*cloud(i,j))*(1.-albedo_i)
#endif

#else
#ifdef NCEPSF
        swr_o(i,j)=dswrf(i,j)*(1.-albedo_w)
#else
        swr_o(i,j)=tot_sol(i,j)*(1.-0.62*cloud(i,j))*(1.-albedo_w)
#endif
#endif

#ifdef ICE
!   conductivity
        if(aice(i,j).gt.0.d0)then
           h(i,j)=volice(i,j)/(areat(i,j,1)*aice(i,j))*1.d-2
           if(h(i,j).gt.4.0d0) then
             h(i,j)=4.d0+(h(i,j)-4.0d0)*26.d0
           endif
           conduc(i,j)=aice(i,j)*k_i*1.d-4
     &        *(t(i,j,1)-t_top(i,j))/h(i,j)
        else
           conduc(i,j) = 0.d0
        endif

!   net flux
        netq_i(i,j)=swr_i(i,j)+lwr_i+senh_i+lath_i+conduc(i,j)

!   change for melting condition
        if(netq_i(i,j).gt.0.d0) then
#ifdef NCEPSF
            swr_i(i,j)=aice(i,j)*dswrf(i,j)*(1.-albedo_melt)
#else
            swr_i(i,j)=aice(i,j)*tot_sol(i,j)
     &       *(1.-0.62*cloud(i,j))*(1.-albedo_melt)
#endif
        endif

        if(aice(i,j).gt.0.d0)then
           netq_i(i,j)=swr_i(i,j)+lwr_i+senh_i+lath_i+conduc(i,j)
        else
           netq_i(i,j)=0.d0
        endif

        netq_o(i,j)=swr_o(i,j)-conduc(i,j)
     &   +(1.-aice(i,j))*(lwr_o(i,j)
     &   +alp_h(i,j)*senh_o(i,j)+alp_w(i,j)*lath_o(i,j))
#else
        netq_o(i,j)=swr_o(i,j)+lwr_o(i,j)+senh_o(i,j)+lath_o(i,j)
#endif

        fricv3(i,j)=(1.d6*(win(i,j)**3)
     &   +6.d0*stdev_wind(i,j)**2*win(i,j)*1.d2)
     &       *(rho_a*cd_a_o/rho_w)**(1.5)

#ifdef DEBUG1
      if(nkai.ge.1670 .and. nkai .le. 1677 .and. ip.eq.14
     & .and. i.eq.16 .and. j.eq.25) then
        write(*,*) nkai,'fricv3',aice(i,j),win(i,j),stdev_wind(i,j),
     &    rho_a,cd_a_o,rho_w
      endif
#endif

#ifdef EVAPO
        evapo=1.d2*rho_a*c_la_w*win(i,j)*(q_a(i,j)-q_w)/rho_w
        evapo=dmin1(evapo,0.d0)
        wflux(i,j)=wflux(i,j)+alp_w(i,j)*evapo
#endif

      endif
      enddo
      enddo

#ifdef ICE
      call ice_pre
#endif
c
#ifdef EVAPO
#ifndef NESTED
! adjust mean sea level
      dum(1) = 0.
      dum(2) = 0.
      do j=3,jml-2
      do i=3,iml-2
        dum(1)=dum(1)+areat(i,j,1)
        dum(2)=dum(2)+hcl(i,j)*areat(i,j,1)
     &       +volice(i,j)*rho_i/rho_w
      enddo
      enddo

      call sum_all(dum,2)
      totarea = dum(1)
      htmnpice = dum(2)
!      if(ip.eq.imaster .and. mod(nkai,nmonth).eq.1) then
!        write(*,*) 'adjust of water flux at ',nkai,
!     &    htmnpice/totarea,htmnpice/totarea*3600,'(cm/day)'
!      endif
! adjustment timescale 1yes
      tscalewf = 1./3600/360
      do j=3,jml-2
      do i=3,iml-2
        wflux(i,j) = wflux(i,j) - htmnpice/totarea*tscalewf
      enddo
      enddo
#endif
#endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#else
c
      real*8 awind,xwind,ywind,atemp,rhmdt,preci,apres,swrad,
     &    fca,fcb,e_w_d,e_w,q_w,lath,senh,lwrad,dq,evapo
      integer nsec,nhourb,nhoura
c
#ifdef COUPLE
      real hnet,srad,evpr,prcp,snow,taux,tauy,tkef
       common /cplflux/ hnet(im,jm),srad(im,jm),evpr(im,jm),
     &  prcp(im,jm),snow(im,jm),taux(im,jm),tauy(im,jm),
     &  radl(im,jm),hsns(im,jm),hltn(im,jm),
     &  tkef(im,jm),ncall1,ncall2
c
      do j = 3,jml-2
      do i = 3,iml-2
        wflux(i,j) = -prcp(i,j)-snow(i,j)+evpr(i,j)
        swr_o(i,j) = srad(i,j)
        netq_o(i,j)= hnet(i,j)
        fricv3(i,j)= tkef(i,j)
        wsx(i,j)   = taux(i,j)
        wsy(i,j)   = tauy(i,j)
      enddo
      enddo
#else
c
ccc   interpolate atmospheric condition
c
      nsec=(nkai-1)*nint(dtuv)
      nhourb=1+nsec/(60*60)
      nhoura=nhourb+1
      fca=dble(nsec-(nhourb-1)*60*60)/(60.*60.)
      fcb=1.-fca
c
      awind=(fcb*awind_h(nhourb)+fca*awind_h(nhoura))*1.d-2
      xwind=(fcb*xwind_h(nhourb)+fca*xwind_h(nhoura))
      ywind=(fcb*ywind_h(nhourb)+fca*ywind_h(nhoura))
      atemp=(fcb*atemp_h(nhourb)+fca*atemp_h(nhoura))
      rhmdt=(fcb*rhmdt_h(nhourb)+fca*rhmdt_h(nhoura))
      preci=(fcb*preci_h(nhourb)+fca*preci_h(nhoura))
      apres=(fcb*apres_h(nhourb)+fca*apres_h(nhoura))
      swrad=(fcb*swrad_h(nhourb)+fca*swrad_h(nhoura))
      uref=(fcb*uref_h(nhourb)+fca*uref_h(nhoura))
      vref=(fcb*vref_h(nhourb)+fca*vref_h(nhoura))
c
c      if(ip.eq.imaster)then
c      write(*,*)nkai
c      write(*,*)fca,fcb,nhourb,nhoura
c      write(*,*)awind,xwind,ywind
c      write(*,*)atemp,rhmdt,preci,apres,swrad
c         write(*,*)uref,vref
c      endif
c      stop
c
ccc   calculate surface forcing
      do j=3,jml-2
      do i=3,iml-2
      if(tex(i,j,1).eq.1.d0) then
c
ccc   sensible heat
        senh=1.d-4*rho_a*c_p_a*c_stab*awind*(atemp-t(i,j,1))
c
ccc   latant heat and evapolation
        e_w=a_w*exp((b_w-t(i,j,1)/d_w)
     &      *t(i,j,1)/(t(i,j,1)+c_w))
        q_w=0.9815*0.622*e_w/(apres-0.378*e_w)
        dq=q_w*(rhmdt-1)*(1-q_w)/(1-q_w+rhmdt*q_w)
        lath=1.d-4*rho_a*l_w*c_la_w*awind*dq
        lath=dmin1(lath,0.d0)
        evapo=1.d2*rho_a*c_la_w*awind*dq/rho_w
c
ccc  longwave radiation
        e_w_d=apres*(dq+q_w)/(0.622+0.378*(dq+q_w))
        lwrad=-1.d-4*( eps_w*cboltz*(atemp+273.15)**4
     &        *(0.39-0.05*sqrt(e_w_d))
     &        +4.*eps_w*cboltz*(atemp+273.15)**3
     &        *(t(i,j,1)-atemp) )
c
ccc   short wave radiation
        swr_o(i,j)=swrad
c
ccc   net heat flux to ocean
        netq_o(i,j)=swrad+lwrad+senh+lath
c
ccc   surface frictional velocity **3
        fricv3(i,j)=1.d6*awind**3
     &       *(rho_a*cd_a_o/rho_w)**(1.5)
c
ccc   water flux to ocean
        wflux(i,j)=preci+evapo
c
        if(i.eq.iw1d.and.ls(ip)+j-1.eq.jw1d)then
           wrw(1)=wflux(i,j)
           wrw(2)=preci
           wrw(3)=evapo
        endif
c
      else
        wflux(i,j)=0.d0
        swr_o(i,j)=0.d0
        netq_o(i,j)=0.d0
        fricv3(i,j)=0.d0
      endif
c
c
      if(ex(i,j,1).eq.1.d0)then
ccc   wind stress
c        wsx(i,j)=rho_a*1.d-3*cd_a_o*(xwind-u(i,j,1))
c     &       *dsqrt((xwind-u(i,j,1))**2+(ywind-v(i,j,1))**2)
c        wsy(i,j)=rho_a*1.d-3*cd_a_o*(ywind-v(i,j,1))
c     &       *dsqrt((xwind-u(i,j,1))**2+(ywind-v(i,j,1))**2)
         wsx(i,j)=rho_a*1.d-3*cd_a_o*xwind
     &        *dsqrt(xwind**2+ywind**2)
c     &        *10.d0
         wsy(i,j)=rho_a*1.d-3*cd_a_o*ywind
     &        *dsqrt(xwind**2+ywind**2)
c     &        *10.d0
      else
         wsx(i,j)=0.d0
         wsy(i,j)=0.d0
      endif
c
c
#ifdef DEBUGB
      if(write_1d.ne.0.and.matsno.eq.0)then
         if(i.eq.ii.and.ls(ip)+j-1.eq.jj)then
            write(*,*)nkai,i,ls(ip)+j-1
            write(*,*)'wsx/wsy: ',wsx(i,j),wsy(i,j)
            write(*,*)'netq/sw: ',netq_o(i,j),swr_o(i,j)
            write(*,*)'wat/u*3: ',wflux(i,j),fricv3(i,j)
            write(*,*)'sh/lh/lw: ',senh,lath,lwrad
         endif
      endif
#endif
c
      enddo
      enddo
c
c
ccc   rain temperature
      raintemp=atemp
c
c
      return
      end
c
#endif
#endif