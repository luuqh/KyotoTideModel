!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   IMPORT RESULT FROM NAOVARJ2 TO ASSIMILATION SYSTEM
!
!   4DVAR output code :  Luu Quang Hung
!   Update            :  2009.07.02
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
! $ main - Luu Quang Hung
!-------------------------------------------------------------------

      MODULE mod_nao

      use param
      implicit none

      private
      public loadtidepara,waveheight,equitide,extracttide

! Global variable
      real*8 a(im,jm,tcm),p(im,jm,tcm),h(im,jm),f(tcm)
      real*8 londeg(im,jm),latdeg(im,jm)
      real*8 h0,s0,p0 !astronomical argument by Schwiderski 1980
      integer i,j

!      iyear  = 2007 ! year
!      imon   =    1 ! month
!      iday   =    1 ! day
!      ihour  =    0 ! hour
!      imin   =    0 ! minute

! Convert to tidal time
!      call mjdymd (time,iyear,imon,iday,ihour,imin,0,1)

! Load variables and store into A,P
!      call loadtidepara (im,jm,tm,a,p)

! Calculate them into H
!      call waveheight (time,im,jm,tm,a,p,h)

! $ vset - Matsumoto
!------------------------------------------------------------------
!      BLOCK data vset
!------------------------------------------------------------------

! * TD - UT values are stored in 'etmut'.
!  They are currently filled up until 2030, but they are not
!  correct beyond 1997. It is preferable to replace them by
!  correct values when they are published.

      real(8) :: eqasmj(16),
     +           eqasmi(33),
     +           eqalmj( 7),
     +           eqalmi( 5),
     +           etmut(130),
     +           fr(3,6),
     +           iagsmj(7,16),
     +           iagsmi(7,33),
     +           iaglmj(7, 7),
     +           iaglmi(7, 5)
      real(8) :: pi,rad,deg
      character*3 wns(16),wnl(7)
      character*4 wnsm(33), wnlm(5)



! Major diurnal

      data ((iagsmj(i,j),i=1,7),j=1,7)
     +         /1,-2, 0, 1, 0, 0,-90, ! Q1       1
     +          1,-1, 0, 0, 0, 0,-90, ! O1       2
     +          1, 0, 0, 1, 0, 0, 90, ! M1       3
     +          1, 1,-2, 0, 0, 0,-90, ! P1       4
     +          1, 1, 0, 0, 0, 0, 90, ! K1       5
     +          1, 2, 0,-1, 0, 0, 90, ! J1       6
     +          1, 3, 0, 0, 0, 0, 90/ ! OO1      7
!               m  s  h  p  n  ps

! Major semi-diurnal

      data ((iagsmj(i,j),i=1,7),j=8,16)
     +         /2,-2, 0, 2, 0, 0,  0, ! 2N2      8
     +          2,-2, 2, 0, 0, 0,  0, ! Mu2      9
     +          2,-1, 0, 1, 0, 0,  0, ! N2      10
     +          2,-1, 2,-1, 0, 0,  0, ! Nu2     11
     +          2, 0, 0, 0, 0, 0,  0, ! M2      12
     +          2, 1, 0,-1, 0, 0,180, ! L2      13
     +          2, 2,-3, 0, 0, 1,  0, ! T2      14
     +          2, 2,-2, 0, 0, 0,  0, ! S2      15
     +          2, 2, 0, 0, 0, 0,  0/ ! K2      16
!               m  s  h  p  n  ps

! Diurnal Minor 17
      data ((iagsmi(i,j),i=1,7),j=1,17)
     +         /1,-3, 0, 2, 0, 0,-90, ! 2Q1     1
     +          1,-3, 2, 0, 0, 0,-90, ! Sigma1  2
     +          1,-2, 0, 1,-1, 0,-90, ! Q1'     3
     +          1,-2, 2,-1, 0, 0,-90, ! Rho1    4
     +          1,-1, 0, 0,-1, 0,-90, ! O1'     5
     +          1,-2, 2, 0, 0, 0, 90, ! Tau1    6
     +          1, 0, 0, 1, 1, 0, 90, ! M1'     7
     +          1, 0, 2,-1, 0, 0, 90, ! Kai1    8
     +          1, 1,-3, 0, 0, 1,-90, ! Pi1     9
     +          1, 1,-2, 0,-1, 0, 90, ! P1'     10
     +          1, 1, 0, 0,-1, 0,-90, ! K1'     11
     +          1, 1, 0, 0, 1, 0, 90, ! K1'     12
     +          1, 1, 1, 0, 0,-1, 90, ! Psi1    13
     +          1, 1, 2, 0, 0, 0, 90, ! Phi1    14
     +          1, 2,-2, 1, 0, 0, 90, ! Theta1  15
     +          1, 2, 0,-1, 1, 0, 90, ! J1'     16
     +          1, 3, 0, 0, 1, 0, 90/ ! OO1'    17
!               m  s  h  p  n  ps

! Semi-diurnal minor 16
      data ((iagsmi(i,j),i=1,7),j=18,33)
     +         /2,-3, 2, 1, 0, 0,  0, ! Eps.2   18
     +          2,-2, 2, 0,-1, 0,180, ! Mu2'    19
     +          2,-1, 0, 1,-1, 0,180, ! N2'     20
     +          2,-1, 2,-1,-1, 0,180, ! Nu2'    21
     +          2, 0,-2, 2, 0, 0,180, ! Gamma2  22
     +          2, 0,-1, 0, 0, 1,180, ! Alpha2  23
     +          2, 0, 0, 0,-1, 0,180, ! M2'     24
     +          2, 0, 1, 0, 0,-1,  0, ! Beta2   25
     +          2, 0, 2, 0, 0, 0,  0, ! Delata2 26
     +          2, 1,-2, 1, 0, 0,180, ! Lambda2 27
     +          2, 2,-2, 0,-1, 0,  0, ! S2'     28
     +          2, 2,-1, 0, 0,-1,180, ! R2      29
     +          2, 2, 0, 0, 1, 0,  0, ! K2'     30
     +          2, 3,-2, 1, 0, 0,  0, ! Zeta2   31
     +          2, 3, 0,-1, 0, 0,  0, ! Eta2    32
     +          2, 3, 0,-1, 1, 0,  0/ ! Eta2'   33
!               m  s  h  p  n  ps

! Major long period

      data ((iaglmj(i,j),i=1,7),j=1,7)
     +         /0, 3, 0,-1, 0, 0,  0, ! Mtm     1
     +          0, 2, 0, 0, 0, 0,  0, ! Mf      2
     +          0, 2,-2, 0, 0, 0,  0, ! MSf     3
     +          0, 1, 0,-1, 0, 0,  0, ! Mm      4
     +          0, 1,-2, 1, 0, 0,  0, ! MSm     5
     +          0, 0, 2, 0, 0, 0,  0, ! Ssa     6
     +          0, 0, 1, 0, 0,-1,  0/ ! Sa      7
!               m  s  h  p  n  ps

! Long-period nodal modulations
      data ((iaglmi(i,j),i=1,7),j=1,5)
     +         /0, 3, 0,-1, 1, 0,  0, ! Mtm'    1
     +          0, 2, 0, 0, 1, 0,  0, ! Mf'     2
     +          0, 1, 0,-1, 1, 0,180, ! Mm'     3
     +          0, 1, 0,-1,-1, 0,180, ! Mm'     4
     +          0, 0, 2, 0, 1, 0,180/ ! Ssa'    5
!               m  s  h  p  n  ps

      data fr/280.4606184d0,  36000.7700536d0,  0.00038793d0,
     +        218.3166560d0, 481267.8813420d0, -0.00133000d0,
     +        280.4664490d0,  36000.7698220d0,  0.00030360d0,
     +         83.3532430d0,   4069.0137110d0, -0.01032400d0,
     +        234.9554440d0,   1934.1361850d0, -0.00207600d0,
     +        282.9373480d0,      1.7195330d0,  0.00045970d0/

      data eqasmj/0.072136d0, 0.376763d0, 0.029631d0, 0.175307d0, ! 1- 4
     +            0.529876d0, 0.029630d0, 0.016212d0, 0.023009d0, ! 5- 8
     +            0.027768d0, 0.173881d0, 0.033027d0, 0.908184d0, ! 9-12
     +            0.025670d0, 0.024701d0, 0.422535d0, 0.114860d0/ !13-16

      data eqasmi/0.009545d0, 0.011520d0, 0.013607d0, 0.013702d0, ! 1- 4
     +            0.071081d0, 0.004914d0, 0.005946d0, 0.005667d0, ! 5- 8
     +            0.010251d0, 0.001973d0, 0.010492d0, 0.071886d0, ! 9-12
     +            0.004145d0, 0.007545d0, 0.005666d0, 0.005875d0, !13-16
     +            0.010385d0, 0.006709d0, 0.001037d0, 0.006484d0, !17-20
     +            0.001232d0, 0.002728d0, 0.003123d0, 0.033885d0, !21-24
     +            0.002749d0, 0.001066d0, 0.006697d0, 0.000946d0, !25-28
     +            0.003536d0, 0.034240d0, 0.001228d0, 0.006422d0, !29-32
     +            0.002799d0/                                     !33

      data eqalmj/0.029926d0, 0.156303d0, 0.013695d0, 0.082569d0, ! 1- 4
     +            0.015791d0, 0.072732d0, 0.011549d0/             ! 5- 7

      data eqalmi/0.012405d0, 0.064805d0, 0.005358d0, 0.005419d0, ! 1- 4
     +            0.001799d0/                                     ! 5

      data (etmut(i),i=1,50)        !1901-1950
     +           / 1.7D0,    2.9D0,    4.0D0,    5.2D0,    6.3D0,
     +             7.5D0,    8.6D0,    9.7D0,   10.9D0,   12.1D0,
     +            13.4D0,   14.6D0,   15.7D0,   16.8D0,   17.8D0,
     +            18.7D0,   19.6D0,   20.4D0,   21.1D0,   21.7D0,
     +            22.2D0,   23.0D0,   22.9D0,   23.4D0,   23.7D0,
     +            23.9D0,   23.9D0,   23.5D0,   23.7D0,   24.0D0,
     +            24.1D0,   24.4D0,   24.2D0,   24.4D0,   24.3D0,
     +            24.2D0,   24.2D0,   24.6D0,   24.4D0,   24.6D0,
     +            25.5D0,   25.5D0,   26.4D0,   26.7D0,   26.9D0,
     +            27.4D0,   28.2D0,   28.5D0,   29.2D0,   29.9D0/
      data (etmut(i),i=51,100)      !1951-2000
     +           /30.0D0,   30.6D0,   31.0D0,   31.3D0,   31.23D0,
     +            31.50D0,  31.92D0,  32.43D0,  32.91D0,  33.37D0,
     +            33.78D0,  34.246D0, 34.735D0, 35.403D0, 36.151D0,
     +            37.000D0, 37.887D0, 38.756D0, 39.708D0, 40.710D0,
     +            41.689D0, 42.825D0, 43.961D0, 44.998D0, 45.986D0,
     +            47.000D0, 48.038D0, 49.105D0, 50.105D0, 50.981D0,
     +            51.817D0, 52.578D0, 53.438D0, 54.090D0, 54.640D0,
     +            55.117D0, 55.585D0, 56.098D0, 56.574D0, 57.227D0,
     +            57.962D0, 58.545D0, 59.590D0, 60.406D0, 61.250D0,
     +            61.980D0, 62.638D0, 63.284D0, 63.664D0, 64.600D0/
      data (etmut(i),i=101,130)     !2001-2030
     +           /65.000D0, 65.000D0, 66.000D0, 67.000D0, 67.000D0,
     +            68.000D0, 69.000D0, 69.000D0, 70.000D0, 71.000D0,
     +            71.000D0, 72.000D0, 73.000D0, 73.000D0, 74.000D0,
     +            74.000D0, 75.000D0, 76.000D0, 76.000D0, 77.000D0,
     +            78.000D0, 78.000D0, 79.000D0, 80.000D0, 80.000D0,
     +            81.000D0, 82.000D0, 82.000D0, 83.000D0, 84.000D0/

      data wns     /'q1 ','o1 ','m1 ','p1 ','k1 ','j1 ',
     +              'oo1','2n2','mu2','n2 ','nu2','m2 ',
     +              'l2 ','t2 ','s2 ','k2 '/
      data wnl     /'mtm','mf ','msf','mm ','msm','ssa','sa '/
      data wnsm    /'2q1 ','sgm1','q1s ','rho1','o1s ','tau1',
     +              'm1s ','kai1','pi1 ','p1s ','k1s ','k1s2',
     +              'psi1','phi1','the1','j1s ','oo1s','eps2',
     +              'mu2s','n2s ','nu2s','gam2','alp2','m2s ',
     +              'bet2','dlt2','lmd2','s2s ','r2  ','k2s ',
     +              'zet2','et2 ','et2s'/
      data wnlm     /'mtms','mfs ','mms1','mms2','ssas'/

      data pi      /3.14159265358979d0/
      data rad     /1.745329251994329d-2/
      data deg     /5.729577951308232d+1/

      contains

! $ waveheight - Luu Quang Hung
!-------------------------------------------------------------------

      SUBROUTINE waveheight (now_time,h)

!-------------------------------------------------------------------

! Produce heights of tidal wave

      use datetime
      use mod_mpi
      implicit none

#include "common.h"

      type(TDate),intent(in) :: now_time
      real*8,intent(out) :: h(im,jm)

      integer i,j,k
      real*8 time,sc !,f
      integer yy,mm,dd,hh,mn,iflag
      iflag = 1

      yy = now_time%year
      mm = now_time%month
      dd = now_time%day
      hh = now_time%hour
      mn = now_time%min
      sc = now_time%sec
      call mjdymd(time,yy,mm,dd,hh,mn,sc,iflag)
      ! call bcast_dble1(time,1)	  

! Set waveheight at every time step
      do i = 1,im
         do j = 1,jm
            h(i,j) = 0.d0
         end do
      end do

      do k = 1,tcm
         if(ut_constituent(k).eq.1)then

            ! k = 5
            ! f = 24.d0/23.9344d0 ! test K1
            ! call frequency (k,f,yy,time)

            do i = 1,im
               do j = 1,jm
                  h(i,j)=h(i,j)+a(i,j,k)*dcos(2*pi*f(k)*time+p(i,j,k))
               end do
            end do
         endif
      end do

      return
      END SUBROUTINE


	  
! $ waveheight - Luu Quang Hung
!-------------------------------------------------------------------

      SUBROUTINE equitide (now_time,zeta)

!-------------------------------------------------------------------

! Produce heights of tidal wave

      use datetime
      use mod_mpi
      implicit none

#include "common.h"

      type(TDate),intent(in) :: now_time
      real*8,intent(out) :: zeta(im,jm)

      integer i,j,k
      real*8 time,sc !,f
      integer yy,mm,dd,hh,mn,iflag,m
      real*8 kamp,chi,eqa ! Schiderski 1980; Matsumoto et al. 2000
      iflag = 1

      yy = now_time%year
      mm = now_time%month
      dd = now_time%day
      hh = now_time%hour
      mn = now_time%min
      sc = now_time%sec
      call mjdymd(time,yy,mm,dd,hh,mn,sc,iflag)
      ! call bcast_dble1(time,1)	  

      call astro(now_time)
	  
! Set waveheight at every time step
      do i = 1,im
         do j = 1,jm
            zeta(i,j) = 0.d0
         end do
      end do	  
	  
      do k = 1,tcm
         if(ut_constituent(k).eq.1)then

            if ((k.ge.1).and.(k.le.16)) then
               m = k
               eqa = eqasmj(m)
            elseif((k.ge.17).and.(k.le.49)) then
               m = k-16
               eqa = eqasmi(m)
            elseif((k.ge.50).and.(k.le.56)) then
               m = k-49
               eqa = eqalmj(m)
            elseif((k.ge.57).and.(k.le.61)) then
               m = k-56
               eqa = eqalmi(m)
            else
               m = 1
               eqa = 0.
               if(ip==imaster) then
                  write(*,*)'out-of-range in mod_nao.F90/equitide'
               endif
            endif    

            kamp = eqa/3.74298732d0*100. ! use NAO parameters (m->cm)
            chi = 0.d0 ! no astronomical argument (as default or single constituent)
            sigma = 2*pi*f(k)

            select case (k) ! use Schwiderski parameters (forced for K1, O1, M2, S2 experiments)
            case (5)  !k1
               kamp = 14.1565 !cm
               chi = h0 + pi/2.
            case (2)  !o1
               kamp = 10.0514 !cm
               chi = h0 - 2*s0 - pi/2.
            case (4)  !p1
               kamp =  4.6843 !cm
               chi = h0 - pi/2.
            case (1)  !q1
               kamp =  1.9256 !cm
               chi = h0 - 3*s0 - p0 - pi/2.
            case (12) !m2
               kamp = 24.2334 !cm
               chi = 2*h0 - 2*s0
            case (15) !s2
               kamp = 11.2841 !cm
               chi = 0.
            case (10)  !n2
               kamp =  4.6398 !cm
               chi = 2*h0 - 3*s0 + p0
            case (16)  !k2
               kamp =  3.0704 !cm
               chi = 2*h0
            case (51)  ! mf
               kamp =  4.1742 !cm
               chi = 2*s0
            case (53)  ! mm
               kamp =  2.2026 !cm
               chi = s0 - p0
            case (55)  ! ssa
               kamp =  1.9446 !cm
               chi = 2*h0
            case default
               kamp = 0 !no values - should be add for other components, or may instead use above NAO
               chi = 0
            end select
			
            do i = 1,im
               do j = 1,jm
                  if(k.le.7) then !major diurnal
                     zeta(i,j) = zeta(i,j) + kamp*dsin(2*latdeg(i,j))
     &                         * dcos(sigma*time+londeg(i,j)+chi)
                  elseif((k.ge.8).and.(k.le.16)) then !major semidiurnal
                     zeta(i,j) = zeta(i,j) + kamp*dsin(latdeg(i,j))**2
     &                         * dcos(sigma*time+2*londeg(i,j)+chi)
                  elseif((k.ge.50).and.(k.le.56)) then ! minor tides and long tides
                     zeta(i,j) = zeta(i,j) + kamp
     &                  *(3*dsin(latdeg(i,j))**2-2)*dcos(sigma*time+chi)
                  else 
                     zeta(i,j) = zeta(i,j) + 0.
                  endif				  
               end do
            end do
         endif
      end do

      return
      END SUBROUTINE
	  


! $ loadtidepara - Luu Quang Hung
!-------------------------------------------------------------------
      SUBROUTINE loadtidepara
!-------------------------------------------------------------------

! Produce heights of tidal wave
      use mod_mpi
      implicit none

#include "common.h"

      character*80  fname,wname
      real*8 ag(img,jmg),pg(img,jmg),long(img,jmg),latg(img,jmg),fg
      real*4 at,pt,dxdegtide,dydegtide,slontide,slattide
      integer i,j,k,m,wtype
 
      do k = 1,tcm

         if ((k.ge.1).and.(k.le.16)) then
            m = k
            wtype = 1
            wname = trim(wns(m))
         elseif((k.ge.17).and.(k.le.49)) then
            m = k-16
            wtype = 2
            wname = trim(wnsm(m))
         elseif((k.ge.50).and.(k.le.56)) then
            m = k-49
            wtype = 3
            wname = trim(wnl(m))
         elseif((k.ge.57).and.(k.le.61)) then
            m = k-56
            wtype = 4
            wname = trim(wnlm(m))
         else
            m = 1
            wtype = 0
            wname = wns(m)
            if(ip==imaster) then
               write(*,*)'out-of-range in mod_nao.F90/loadtidepara'
            endif
         endif                
         call frequency (m,fg,wtype)
         f(k) = fg

         if(ip==imaster) then
            fname = trim(datadir)//'tide/'//trim(wname)//'.dat'
            write(*,*)'loadtide:',trim(fname),'.'
            open(80,file=trim(fname),form='unformatted')

            do i = 1,img
               do j = 1,jmg
                  read(80)at,pt
                  ag(i,j) = dble(at)
                  pg(i,j) = dble(pt)
               end do
            end do
			
            if(k.eq.51) then ! Mf correction
               do j = 1,jmg
                  do i = 1,int(img/3)
                     ag(i,j) = ag(i,j)*3.2
                     ! pg(i,j) = pg(i,j)*0.8
                  enddo
                  do i = int(img/3)+1,img
                     ag(i,j) = ag(i,j)*0.3
                     ! pg(i,j) = pg(i,j)*2.6
                  enddo
               enddo			
            endif			   
			
            close(80)

         endif !imaster
		 
         call bcast_dble(ag,img*jmg)
         do j = 1,jml
            do i = 1,iml
               a(i,j,k) = ag(lw(ip_x)-1+i,ls(ip_y)-1+j)
            enddo
         enddo

         call bcast_dble(pg,img*jmg)
         do j = 1,jml
            do i = 1,iml
               p(i,j,k) = pg(lw(ip_x)-1+i,ls(ip_y)-1+j)
            enddo
         enddo

      enddo

      if(ip==imaster) then
         dxdegtide = 1./54.
         dydegtide = 1./72.
         slontide = 127.+1./6.+206./18.-3./54.
         slattide = 30.5+192./24.-3./72.
         do i = 1,img
            do j = 1,jmg
               long(i,j) = dble(slontide+dxdegtide*(i-1))/180.d0*pi
               latg(i,j) = dble(slattide+dydegtide*(j-1))/180.d0*pi
            end do
         end do
      endif	  

      call bcast_dble(long,img*jmg)
      do j = 1,jml
         do i = 1,iml
            londeg(i,j) = long(lw(ip_x)-1+i,ls(ip_y)-1+j)
         enddo
      enddo
	  
      call bcast_dble(latg,img*jmg)
      do j = 1,jml
         do i = 1,iml
            latdeg(i,j) = latg(lw(ip_x)-1+i,ls(ip_y)-1+j)
         enddo
      enddo

      return
      END SUBROUTINE loadtidepara



! $ astro - Luu Quang Hung
!-------------------------------------------------------------------
      SUBROUTINE astro(now_time)
!-------------------------------------------------------------------

      use datetime
      type(TDate),intent(in) :: now_time

      integer monthday(12),yy,mm,dd
      real*8 T,D
      data monthday /31,28,31,30,31,30,31,31,30,31,30,31/ !days in a month (not leap year)

      ! astronomical argument by Schwiderski 1980
      yy = now_time%year ! year number
      mm = now_time%month
      if (yy.lt.1975) yy=2000
      dd = 0
      if (mm.gt.1) then
         do i=1,mm
            dd = dd+monthday(i)
         enddo
      endif
      dd = dd +now_time%day
      if (mod(yy,4).eq.0) dd = dd +1 !leap year
      D = dd +365*(yy-1975)+int(365*(yy-1975))
      T = (27392.500528 +1.0000000365*D)/36525
      h0 = 279.696680 +36000.768930485*T +0.000303*T*T
      s0 = 270.434358 +481267.88314137*T -0.001133*T*T 
     &     +1.9*T*T*T*1d-6
      p0 = 334.329653 +4069.0340329575*T -0.010325*T*T 
     &     -1.2*T*T*T*1d-5

      return
      END SUBROUTINE astro



! $ frequency - Luu Quang Hung
!------------------------------------------------------------------
      SUBROUTINE frequency (iw,freq,wtype)
!------------------------------------------------------------------

      implicit none

      real*8,intent(out) :: freq
      integer,intent(in) :: iw,wtype
      real*8,parameter :: r36525 = 1.d0/36525.d0
      real*8 xarg(7,33),ftau,temp
       
      do i = 1,7
         if(wtype.eq.1)then
            do j = 1,16
               xarg(i,j) = iagsmj(i,j)
            enddo
         elseif(wtype.eq.2)then
            do j = 1,33
               xarg(i,j) = iagsmi(i,j)
            enddo
         elseif(wtype.eq.3)then
            do j = 1,7
               xarg(i,j) = iaglmj(i,j)
            enddo
         else
            do j = 1,5
               xarg(i,j) = iaglmi(i,j)
            enddo
         endif
      enddo

      ftau = 360.d0/r36525 - fr(2,2) + fr(2,3)
      temp = ftau*xarg(1,iw)    + fr(2,2)*xarg(2,iw)
     +        + fr(2,3)*xarg(3,iw) + fr(2,4)*xarg(4,iw)
     +        + fr(2,5)*xarg(5,iw) + fr(2,6)*xarg(6,iw)
      temp = temp*r36525
      freq = temp/360.     ! cycles/day

      return
      END subroutine frequency





! $ extractide - Luu Quang Hung
!-------------------------------------------------------------------
      SUBROUTINE extracttide(ax,px,fx)
!-------------------------------------------------------------------

! Produce heights of tidal wave
      use mod_mpi
      implicit none

      real*8,intent(out):: ax(im,jm,tcm),px(im,jm,tcm),fx(tcm)
      integer:: i,j,k

      do k=1,tcm
         fx(k) = f(k)
         do i=1,iml
         do j=1,jml
            ax(i,j,k) = a(i,j,k)
            px(i,j,k) = p(i,j,k)
         enddo
         enddo
      enddo

      return
      END SUBROUTINE extracttide



! $ mjdymd - Matsumoto
!---------------------------------------------------------------------
      SUBROUTINE mjdymd(xmjd  , iy    , im    , id    , ih    ,
     +                  imin  , sec  , iflag                  )
!---------------------------------------------------------------------

! xmjd  : modified julian date
! iy    : year
! im    : month
! id    : day
! ih    : hour
! imin  : minute
! sec  : second
! iflag : 1 -> YMDHMS to MJD
!         2 -> MJD to YMDHMS
! Date must be within the years Mar. 1, 1900 to Feb. 28, 2100

      implicit none

      real*8,parameter :: xjd0 = 2400000.5d0
      real*8,parameter :: half =       0.5d0

      real*8,intent(inout) :: xmjd,sec
      integer,intent(inout) ::  iy,im,id,ih,imin
      integer,intent(in) :: iflag

      real*8 y,xjd,fsec,frc
      integer m,c,nd,e,nf,ifr,mjd,isec


! -----< YMDHMS to MJD >-----

      if (iflag.eq.1) then

         y = dfloat(iy - 1)

         if (im.gt.2) then
            m = im
            y = y + 1
         else
            m = im + 12
         endif

         xjd  = int(365.25d0*y) + int(30.6001d0*(m+1)) - 15
     +        + 1720996.5d0     + id
         xmjd = xjd - xjd0

         fsec = dfloat(ih)*3600.d0 + dfloat(imin)*60.d0 + sec

         xmjd = xmjd + fsec/86400.d0

! -----< MJD to YMDHMS >-----

      elseif (iflag.eq.2) then

         mjd  = xmjd
         xjd  = dfloat(mjd) + xjd0
         c    = int(xjd + half) + 1537
         nd   = int((c - 122.1d0)/365.25d0 )
         e    = int(365.25d0*nd)
         nf   = int((c - e)/30.6001d0)

         ifr  = int(xjd + half)
         frc  = xjd + half - dfloat(ifr)
         id   = c - e - int(30.6001d0*nf) + frc
         im   = nf - 1 - 12*int(nf/14)
         iy   = nd - 4715 - int((7+im)/10)

         sec  = (xmjd-dfloat(mjd))*86400.d0
         isec = sec
         if ((sec-isec).gt.0.5d0) isec = isec + 1
         ih   = isec/3600
         imin = (isec - ih*3600)/60
         isec = isec - ih*3600 - imin*60

      else

         print*,'!!! Error in <mjdymd>. iflag should be 1 or 2.'
         stop

      endif

      return
      END subroutine

      END MODULE mod_nao
c ----------------------< End of program >----------------------
c \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/


