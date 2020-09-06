!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     subroutine brtrop                                             c

!  $Id: brtrop.F90 13 2008-12-12 10:53:49Z ishikawa $

!     In this subroutine, the barotropic velocity component (depth- c
!     avergaged) is obtained by solving shallow-water equation,     c
!     using time-split method.                                      c
!                                                                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine brtrop

      use param
      use mod_mpi
#ifdef BNDTIDE
      use mod_nao
      use mod_dump
      use mod_time
#endif
      implicit none
#include "common.h"

      real*8 uvis,vvis,viscoe,flmd,fphi,
     &     uvis_l,vvis_l,uvis_b,vvis_b,fclfr,uvlimit
#ifdef NESTED
      real*8 fbcoef,factu
#endif
      integer matsu1,matsu2,i,j,k,n,jsn,jen
      real*8 umm(im,jm),vmm(im,jm)
#ifdef DEBUG11X
      integer,save :: ncnt_bt=0
      integer :: ncstat
#endif
#ifdef BNDTIDE
      real*8 h_tide(im,jm),zeta(im,jm)
      real*8 vel,gae,dump_l,dump_r,dump_s,dump_h,eta_drag,frac
      type(Tdate) :: date_trop
      ! real*8,parameter :: f1 = 1./36.d4/10. 
      ! real*8,parameter :: f2 = 1.d6 !cm-1
      real*8,parameter :: alpha = 0.94 !0.9 !0.94
      real*8,parameter :: beta = 0.736 !0.69 !0.736
      real*8,parameter :: f1 = 5.5d-6
      real*8,parameter :: f1h = f1*10.
#endif

!    facters for calm spin-up
!    fclfr : laplacian friction
#ifdef BISMFRIC
      fclfr=0.d0
#else
      fclfr=1.d0
#endif
#ifdef PC68M
!      fclfr=dmax1( dble(nxmonth-nkai)/dble(nxmonth)*10.d0,0.d0 )
!      fclfr=dmax1( dble(nxmonth-nkai)/dble(nxmonth)*10.d0,1.d0 )
      fclfr=1.d0
#endif

#ifdef NWNPAC
      fclfr=0.d0
#endif
#ifdef MASKLVIS
      fclfr=1.d0
#endif
#ifdef BNDTIDE
       fclfr = 1.d0
#endif


      do j = 1,jml
      do i = 1,iml
!   2nd of mastuno step in subroutine tracli
        if(matsn2.eq.1) then
           um(i,j)=umd(i,j)
            vm(i,j)=vmd(i,j)
            sfun(i,j) = sfund(i,j)
            sfvn(i,j) = sfvnd(i,j)
!   ordinary time step or 1st of matsuno in subroutine tracli
        else
            umd(i,j)=um(i,j)
            vmd(i,j)=vm(i,j)
            sfund(i,j) = sfun(i,j)
            sfvnd(i,j) = sfvn(i,j)
        endif

!   ht is always start from hclb
         ht(i,j) = hclb(i,j)

!   initialize
         umstar(i,j)=0.d0
         vmstar(i,j)=0.d0

         htb(i,j) = ht(i,j)
         umb(i,j) = um(i,j)
         vmb(i,j) = vm(i,j)
         sfunb(i,j) = sfun(i,j)
         sfvnb(i,j) = sfvn(i,j)

         umm(i,j) = um(i,j)/dble(ntrop*2+1)
         vmm(i,j) = vm(i,j)/dble(ntrop*2+1)
      enddo
      enddo

#ifdef BNDTIDE
      date_trop = date_now%tot_sec	  
      call equitide(date_trop,zeta)	  
#endif


#ifdef NESTED

#ifdef REST_INI
#ifdef ROKKA
      factu=dmin1(dble(nkai)/dble(nxmonth),1.d0)
#endif
#ifdef NWNPAC
      factu=dmin1(dble(nkai)/dble(nxmonth),1.d0)
#endif
#ifdef JP68M
      factu=dmin1(dble(nkai)/dble(nxmonth),1.d0)
#endif
#else
      factu=1.d0
#endif

! west
      if(ip_x .eq. 0) then
      do j=1,jml
         do i=1,indg+2
            htwbt(i,j)=hclbwbt(i,j)
            htbwbt(i,j)=hclbwbt(i,j)
            wtawbt(i,j)=(wtwbc(i,j,mb1)
     $           +(wtwbc(i,j,mb2)-wtwbc(i,j,mb1))*cbf)*factu
            wspwbt(i,j)=wtbwbt(i,j)
         enddo
         do i=1,2
            umwbt(i,j,1)=umb(i+2,j)
            umwbt(i,j,2)=(umwbc(i,j,mb1)
     $           +(umwbc(i,j,mb2)-umwbc(i,j,mb1))*cbf)*factu
            vmwbt(i,j,1)=vmb(i+2,j)
            vmwbt(i,j,2)=(vmwbc(i,j,mb1)
     $           +(vmwbc(i,j,mb2)-vmwbc(i,j,mb1))*cbf)*factu
         enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do j=1,jml
         do i=1,indg+2
            htebt(i,j)=hclbebt(i,j)
            htbebt(i,j)=hclbebt(i,j)
            wtaebt(i,j)=(wtebc(i,j,mb1)
     $           +(wtebc(i,j,mb2)-wtebc(i,j,mb1))*cbf)*factu
            wspebt(i,j)=wtbebt(i,j)
         enddo
         do i=1,2
            umebt(i,j,1)=umb(iml-2-i,j)
            umebt(i,j,2)=(umebc(i,j,mb1)
     $           +(umebc(i,j,mb2)-umebc(i,j,mb1))*cbf)*factu
            vmebt(i,j,1)=vmb(iml-2-i,j)
            vmebt(i,j,2)=(vmebc(i,j,mb1)
     $           +(vmebc(i,j,mb2)-vmebc(i,j,mb1))*cbf)*factu
         enddo
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do i=1,iml
         do j=1,jndg+2
            htsbt(i,j)=hclbsbt(i,j)
            htbsbt(i,j)=hclbsbt(i,j)
            wtasbt(i,j)=(wtsbc(i,j,mb1)
     $           +(wtsbc(i,j,mb2)-wtsbc(i,j,mb1))*cbf)*factu
            wspsbt(i,j)=wtbsbt(i,j)
         enddo
         do j=1,2
            umsbt(i,j,1)=umb(i,j+2)
            vmsbt(i,j,1)=vmb(i,j+2)
            umsbt(i,j,2)=(umsbc(i,j,mb1)
     $           +(umsbc(i,j,mb2)-umsbc(i,j,mb1))*cbf)*factu
            vmsbt(i,j,2)=(vmsbc(i,j,mb1)
     $           +(vmsbc(i,j,mb2)-vmsbc(i,j,mb1))*cbf)*factu
         enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do i=1,iml
         do j=1,jndg+2
            htnbt(i,j)=hclbnbt(i,j)
            htbnbt(i,j)=hclbnbt(i,j)
            wtanbt(i,j)=(wtnbc(i,j,mb1)
     $           +(wtnbc(i,j,mb2)-wtnbc(i,j,mb1))*cbf)*factu
            wspnbt(i,j)=wtbnbt(i,j)
         enddo
         do j=1,2
            umnbt(i,j,1)=umb(i,jml-2-j)
            vmnbt(i,j,1)=vmb(i,jml-2-j)
            umnbt(i,j,2)=(umnbc(i,j,mb1)
     $           +(umnbc(i,j,mb2)-umnbc(i,j,mb1))*cbf)*factu
            vmnbt(i,j,2)=(vmnbc(i,j,mb1)
     $           +(vmnbc(i,j,mb2)-vmnbc(i,j,mb1))*cbf)*factu
         enddo
      enddo
      endif
#endif

	
!   loop for time-spliting scheme

      do 1200 n=1,ntrop*2

#ifdef NESTED
      fbcoef=dble(n)/dble(ntrop)
#endif

      if(mod(n,modmtnb).eq.1) then
         c2dttr = dttr
         matsu1 = 1
         matsu2 = 1
c
         do j = 1,jml
         do i = 1,iml
            umb(i,j) = um(i,j)
            vmb(i,j) = vm(i,j)
            htb(i,j) = ht(i,j)
            sfunb(i,j) = sfun(i,j)
            sfvnb(i,j) = sfvn(i,j)
         enddo
         enddo
#ifdef NESTED
! west
      if(ip_x .eq. 0) then
         do j=1,jml
         do i=1,indg+2
            htbwbt(i,j)=htwbt(i,j)
         enddo
         enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
         do j=1,jml
         do i=1,indg+2
            htbebt(i,j)=htebt(i,j)
         enddo
         enddo
      endif
! south
      if(ip_y .eq. 0) then
         do i=1,iml
         do j=1,jndg+2
            htbsbt(i,j)=htsbt(i,j)
         enddo
         enddo
       endif
! north
      if(ip_y .eq. jpe-1) then
         do i=1,iml
         do j=1,jndg+2
            htbnbt(i,j)=htnbt(i,j)
         enddo
         enddo
      endif
#endif

      else
         c2dttr = dttr2
         matsu1 = 0
         matsu2 = 0
      endif

c
#ifdef BISMFRIC
      call exch_t2d_n2p(sfunb,sfvnb,6)
      call exch_t2d_s2p(sfunb,sfvnb,7)
      call wait_t2d_n2p(sfunb,sfvnb,6)
      call wait_t2d_s2p(sfunb,sfvnb,7)
      call exch_t2d_e2p(sfunb,sfvnb,8)
      call exch_t2d_w2p(sfunb,sfvnb,9)
#endif
c
      do j=1,jml
         dtcor(j)=.5d0*cor(j)*c2dttr
         dtcor2(j)=1.d0+dtcor(j)*dtcor(j)
      enddo
c
#ifdef BISMFRIC
      call wait_t2d_e2p(sfunb,sfvnb,8)
      call wait_t2d_w2p(sfunb,sfvnb,9)
#endif


!   return point for matsuno scheme
 1300  continue

#ifdef BISMFRIC
!   biharmonic Smagorinsky friction scheme
      do j = 2,jml-1
      do i = 2,iml-1
         d2sud2x(i,j)=dydxr*csr(j)
     &        *ex(i-1,j,1)*ex(i,j,1)*ex(i+1,j,1)*(
     &        bsmage_0(i,j)*(sfunb(i+1,j)-sfunb(i,j))
     &        -bsmagw_0(i,j)*(sfunb(i,j)-sfunb(i-1,j)))*areaur(j)
         d2sud2y(i,j)=dxdyr*cs(j)
     &        *ex(i,j-1,1)*ex(i,j,1)*ex(i,j+1,1)*(
     &        bsmagn_0(i,j)*(sfunb(i,j+1)-sfunb(i,j))
     &        -bsmags_0(i,j)*(sfunb(i,j)-sfunb(i,j-1)))*areaur(j)
         d2svd2x(i,j)=dydxr*csr(j)
     &        *ex(i-1,j,1)*ex(i,j,1)*ex(i+1,j,1)*(
     &        bsmage_0(i,j)*(sfvnb(i+1,j)-sfvnb(i,j))
     &        -bsmagw_0(i,j)*(sfvnb(i,j)-sfvnb(i-1,j)))*areaur(j)
         d2svd2y(i,j)=dxdyr*cs(j)
     &        *ex(i,j-1,1)*ex(i,j,1)*ex(i,j+1,1)*(
     &        bsmagn_0(i,j)*(sfvnb(i,j+1)-sfvnb(i,j))
     &        -bsmags_0(i,j)*(sfvnb(i,j)-sfvnb(i,j-1)))*areaur(j)
      enddo
      enddo
#endif

      do j=3,jml-2
      do i=3,iml-2
#ifdef BISMFRIC
!   biharmonic smagorinsky friction scheme
         uvis_b=(-dydxr*csr(j)*(bsmage_0(i,j)
     &        *hrumin(i+1,j)*(d2sud2x(i+1,j)-d2sud2x(i,j))
     &        -bsmagw_0(i,j)
     &        *hrumin(i,j)*(d2sud2x(i,j)-d2sud2x(i-1,j)))
     &        -dxdyr*cs(j)*( bsmagn_0(i,j)
     &        *hrvmin(i,j+1)*(d2sud2y(i,j+1)-d2sud2y(i,j))
     &        -bsmags_0(i,j)
     &        *hrvmin(i,j)*(d2sud2y(i,j)-d2sud2y(i,j-1)))
     &        )*areaur(j)*ex(i,j,1)
         vvis_b=(-dydxr*csr(j)*( bsmage_0(i,j)
     &        *hrumin(i+1,j)*(d2svd2x(i+1,j)-d2svd2x(i,j))
     &        -bsmagw_0(i,j)
     &        *hrumin(i,j)*(d2svd2x(i,j)-d2svd2x(i-1,j)))
     &        -dxdyr*cs(j)*( bsmagn_0(i,j)
     &        *hrvmin(i,j+1)*(d2svd2y(i,j+1)-d2svd2y(i,j))
     &        -bsmags_0(i,j)
     &        *hrvmin(i,j)*(d2svd2y(i,j)-d2svd2y(i,j-1)))
     &        )*areaur(j)*ex(i,j,1)
#else
         uvis_b=0.
         vvis_b=0.
#endif

!   laplacian friction scheme
         viscoe=2.d0*dydxr*(aindm1(i,j)*dhru(i,j)
     &        +(1.-aindm1(i+1,j))*dhru(i+1,j))*csr(j)
     &        +2.d0*dxdyr*(aindm2(i,j)*dhrv(i,j)*cst(j)
     &        +(1.-aindm2(i,j+1))*dhrv(i,j+1)*cst(j+1))
#ifdef MASKLVIS
         uvis_l=fclfr*hduv*(
     &        dydxr*csr(j)*(
     $        5.d-1*(mask_lvis(i-1,j)+mask_lvis(i,j))
     $        *hrumin(i,j)*(sfunb(i-1,j)-sfunb(i,j))
     $        -5.d-1*(mask_lvis(i,j)+mask_lvis(i+1,j))
     &        *hrumin(i+1,j)*(sfunb(i,j)-sfunb(i+1,j)))
     &        +dxdyr*(
     $        5.d-1*(mask_lvis(i,j-1)+mask_lvis(i,j))
     $        *hrvmin(i,j)*cst(j)*(sfunb(i,j-1)-sfunb(i,j))
     &        -5.d-1*(mask_lvis(i,j)+mask_lvis(i,j+1))
     $        *hrvmin(i,j+1)*cst(j+1)*(sfunb(i,j)-sfunb(i,j+1)))
     &        -viscoe*mask_lvis(i,j)*sfunb(i,j) )*areaur(j)

         vvis_l=fclfr*hduv*(
     $        5.d-1*(mask_lvis(i-1,j)+mask_lvis(i,j))
     &        *dydxr*csr(j)*(hrumin(i,j)*(sfvnb(i-1,j)-sfvnb(i,j))
     $        -5.d-1*(mask_lvis(i,j)+mask_lvis(i+1,j))
     &        *hrumin(i+1,j)*(sfvnb(i,j)-sfvnb(i+1,j)))
     &        +dxdyr*(
     $        5.d-1*(mask_lvis(i,j-1)+mask_lvis(i,j))
     $        *hrvmin(i,j)*cst(j)*(sfvnb(i,j-1)-sfvnb(i,j))
     &        -5.d-1*(mask_lvis(i,j)+mask_lvis(i,j+1))
     &        *hrvmin(i,j+1)*cst(j+1)*(sfvnb(i,j)-sfvnb(i,j+1)))
     &        -viscoe*mask_lvis(i,j)*sfvnb(i,j) )*areaur(j)
#else
         uvis_l=fclfr*hduv*(
     &        dydxr*csr(j)*(hrumin(i,j)*(sfunb(i-1,j)-sfunb(i,j))
     &        -hrumin(i+1,j)*(sfunb(i,j)-sfunb(i+1,j)))
     &        +dxdyr*(hrvmin(i,j)*cst(j)*(sfunb(i,j-1)-sfunb(i,j))
     &        -hrvmin(i,j+1)*cst(j+1)*(sfunb(i,j)-sfunb(i,j+1)))
     &        -viscoe*sfunb(i,j) )*areaur(j)

         vvis_l=fclfr*hduv*(
     &        dydxr*csr(j)*(hrumin(i,j)*(sfvnb(i-1,j)-sfvnb(i,j))
     &        -hrumin(i+1,j)*(sfvnb(i,j)-sfvnb(i+1,j)))
     &        +dxdyr*(hrvmin(i,j)*cst(j)*(sfvnb(i,j-1)-sfvnb(i,j))
     &        -hrvmin(i,j+1)*cst(j+1)*(sfvnb(i,j)-sfvnb(i,j+1)))
     &        -viscoe*sfvnb(i,j) )*areaur(j)
#endif

!cc   coriolis semi-implicit

#ifndef EXPTEST
		  
         flmd = cor(j)*vmb(i,j)
     &        +zu(i,j)+uvis_l+uvis_b-grav*dx2r*csr(j)*hrr(i,j)
     &        *(ht(i+1,j+1)+ht(i+1,j)-ht(i,j+1)-ht(i,j))
         fphi = -cor(j)*umb(i,j)
     &        +zv(i,j)+vvis_l+vvis_b-grav*dy2r*hrr(i,j)
     &        *(ht(i+1,j+1)+ht(i,j+1)-ht(i+1,j)-ht(i,j))
 
#else
         flmd = cor(j)*vmb(i,j)
     &        +zu(i,j)+uvis_l+uvis_b-grav*dx2r*csr(j)*hrr(i,j)
     &        *(alpha*(ht(i+1,j+1)+ht(i+1,j)-ht(i,j+1)-ht(i,j))
     &        -beta*(zeta(i+1,j+1)+zeta(i+1,j)-zeta(i,j+1)-zeta(i,j)))
         fphi = -cor(j)*umb(i,j)
     &        +zv(i,j)+vvis_l+vvis_b-grav*dy2r*hrr(i,j)
     &        *(alpha*(ht(i+1,j+1)+ht(i,j+1)-ht(i+1,j)-ht(i,j))
     &        -beta*(zeta(i+1,j+1)+zeta(i,j+1)-zeta(i+1,j)-zeta(i,j)))
#endif


#ifdef BNDTIDE
        dump_l = dump1(i,j)
        dump_r = dump2(i,j)
        dump_s = dump3(i,j)

        vel = sqrt(sfun(i,j)*sfun(i,j)+sfvn(i,j)*sfvn(i,j))
        flmd = flmd - (dump_l*vel+dump_r)*um(i,j)
        fphi = fphi - (dump_l*vel+dump_r)*vm(i,j)

        vel = abs((hrr(i+1,j)-hrr(i-1,j))/dx*sfun(i,j)
     &        +(hrr(i,j+1)-hrr(i,j-1))/dy*sfvn(i,j))/2
        dump_s = dump_s*vel/(hrr(i,j)+(1.d0-ex(i,j,1)))*ex(i,j,1)
        flmd = flmd - dump_s*um(i,j)
        fphi = fphi - dump_s*vm(i,j)

        ! flmd = flmd - cdbar(i,j)*vel*um(i,j) !Schiller and Fiedler 2007
        ! fphi = fphi - cdbar(i,j)*vel*vm(i,j)
		
!        if(lw(ip_x)+i-1 == 164 .and. ls(ip_y)+j-1==7) then
!          write(*,*) um(i,j),umb(i,j),flmd,vel,zu(i,j),uvis_l,uvis_b
!        endif
#endif

         uma(i,j) = ex(i,j,1)*(umb(i,j)+c2dttr/dtcor2(j)
     &        *(flmd+dtcor(j)*fphi))
         vma(i,j) = ex(i,j,1)*(vmb(i,j)+c2dttr/dtcor2(j)
     &        *(fphi-dtcor(j)*flmd))

!   without semi-implicit
!
!          uma(i,j)=ex(i,j,1)*(umb(i,j)+c2dttr*(cor(j)*vm(i,j)
!     $      +zu(i,j)+uvis-grav*dx2r*csr(j)*hrr(i,j)
!     $       *(ht(i+1,j+1)+ht(i+1,j)-ht(i,j+1)-ht(i,j)) ))
!          vma(i,j)=ex(i,j,1)*(vmb(i,j)+c2dttr*(-cor(j)*um(i,j)
!     $      +zv(i,j)+vvis-grav*dy2r*hrr(i,j)
!     $       *(ht(i+1,j+1)+ht(i,j+1)-ht(i+1,j)-ht(i,j)) ))

      enddo
      enddo

#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do j = 1,jml
         do i=1,2
            uma(i+2,j)=(umwbt(i,j,1)+(umwbt(i,j,2)-umwbt(i,j,1))*fbcoef
     &           )*ex(i+2,j,1)
            vma(i+2,j)=(vmwbt(i,j,1)+(vmwbt(i,j,2)-vmwbt(i,j,1))*fbcoef
     &           )*ex(i+2,j,1)
#ifdef BNDTIDE
         uma(i+2,j) = uma(5,j)
         vma(i+2,j) = vma(5,j)
#endif
         enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do j = 1,jml
         do i=1,2
        uma(iml-2-i,j)=(umebt(i,j,1)+(umebt(i,j,2)-umebt(i,j,1))*fbcoef
     &        )*ex(iml-2-i,j,1)
        vma(iml-2-i,j)=(vmebt(i,j,1)+(vmebt(i,j,2)-vmebt(i,j,1))*fbcoef
     &        )*ex(iml-2-i,j,1)
#ifdef BNDTIDE
         uma(iml-2-i,j) = uma(iml-5,j)
         vma(iml-2-i,j) = vma(iml-5,j)
#endif
         enddo
         uma(iml-2,j) = 0.
         vma(iml-2,j) = 0.
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do i=1,iml
         do j=1,2
            uma(i,j+2)=(umsbt(i,j,1)+(umsbt(i,j,2)-umsbt(i,j,1))*fbcoef
     &           )*ex(i,j+2,1)
            vma(i,j+2)=(vmsbt(i,j,1)+(vmsbt(i,j,2)-vmsbt(i,j,1))*fbcoef
     &           )*ex(i,j+2,1)
#ifdef BNDTIDE
         uma(i,j+2) = uma(i,5)
         vma(i,j+2) = vma(i,5)
#endif
         enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do i=1,iml
         do j=1,2
        uma(i,jml-2-j)=(umnbt(i,j,1)+(umnbt(i,j,2)-umnbt(i,j,1))*fbcoef
     &           )*ex(i,jml-2-j,1)
        vma(i,jml-2-j)=(vmnbt(i,j,1)+(vmnbt(i,j,2)-vmnbt(i,j,1))*fbcoef
     &           )*ex(i,jml-2-j,1)
#ifdef BNDTIDE
         uma(i,jml-2-j) = uma(i,jml-5)
         vma(i,jml-2-j) = vma(i,jml-5)
#endif
         enddo
         uma(i,jml-2) = 0.
         vma(i,jml-2) = 0.
      enddo
      endif
#endif

      call exch_t2d_n1p(uma,vma,1)
      call exch_t2d_s1p(uma,vma,2)
      call wait_t2d_n1p(uma,vma,1)
      call wait_t2d_s1p(uma,vma,2)
      call exch_t2d_e1p(uma,vma,10)
      call exch_t2d_w1p(uma,vma,11)

	  
      do j = 3,jml-1
!  um/vm need 1,jml-1 !elevation
      do i = 3,iml-1
         hta(i,j)=htb(i,j)+c2dttr*tex(i,j,1)*.25d0*(
     &        dy*( (ex(i-1,j,1)+ex(i-1,j-1,1))
     &        *((2.d0-ex(i-1,j-1,1))*um(i-1,j)
     &        +(2.d0-ex(i-1,j,1))*um(i-1,j-1))
     &        -(ex(i,j,1)+ex(i,j-1,1))
     &        *((2.d0-ex(i,j-1,1))*um(i,j)
     &        +(2.d0-ex(i,j,1))*um(i,j-1)) )
     &        +dx*( (ex(i,j-1,1)+ex(i-1,j-1,1))
     &        *((2.d0-ex(i-1,j-1,1))*vm(i,j-1)
     &        +(2.d0-ex(i,j-1,1))*vm(i-1,j-1))*cs(j-1)
     &        -(ex(i,j,1)+ex(i-1,j,1))
     &        *((2.d0-ex(i-1,j,1))*vm(i,j)
     &        +(2.d0-ex(i,j,1))*vm(i-1,j))*cs(j) )
     &        )/(areat(i,j,1)+(1.d0-tex(i,j,1)))

#ifdef BNDTIDE
       ! if( abs(ht(i,j)) > 50.) then

       if( hrr(i,j) < 100.d2) then ! if depth is small than 1m then dump in Ishikawa
         ! if (hrr(i,j).gt.0.d0) fdump_h = f2/hrr(i,j)
         ! fdump_h = f1h
         ! hta(i,j) = hta(i,j) - fdump_h*ht(i,j)
       endif

       ! cdbar(i,j)=0. ! (Schiller and Fiedler 2007)
       !  eta_drag=abs(hta(i,j)-htb(i,j))/dttr
       ! frac=.95d-2
       ! if(eta_drag.gt.1d-5)cdbar(i,j)=frac*1.d-4
       ! if((eta_drag.le.1d-5).and.(eta_drag.gt.1d-5)) 
   ! &   cdbar(i,j)=frac*eta_drag*dttr*10.
       ! if(eta_drag.le.1d-7) cdbar(i,j)=frac*1.d-6
#endif

      enddo
      enddo

#ifdef NESTED

#ifdef BNDTIDE

      date_trop = date_trop+diff_dttr
      call waveheight(date_trop,h_tide)
      call equitide(date_trop,zeta)
	  
      call exch_t2d_s1(h_tide,20)
      call wait_t2d_s1(h_tide,20)
      call exch_t2d_w1(h_tide,21)
      call wait_t2d_w1(h_tide,21)
      call exch_t2d_n1(h_tide,22)
      call wait_t2d_n1(h_tide,22)
      call exch_t2d_e1(h_tide,23)
      call wait_t2d_e1(h_tide,23)

#endif !BNDTIDE

! south
      if(ip_y .eq. 0) then
      do i=1,iml
      do j=1,indg+2
         htasbt(i,j)=htbsbt(i,j)+tex(i,j+2,1)*c2dttr*wspsbt(i,j)
      enddo
      enddo
      do i=3,iml-2
      do j=3,jndg+2
#ifdef BNDTIDE
         htsbt(i,j) = h_tide(i,j+2)
#endif
         hta(i,j+2)=(hta(i,j+2)+chfht*c2dttr*(htsbt(i,j)-ht(i,j+2))
     &        )*tex(i,j+2,1)
      enddo
      enddo
      do i=1,iml
      do j=1,2
#ifdef BNDTIDE
         htasbt(i,j) = h_tide(i,j+2)
#endif
         hta(i,j+2)=htasbt(i,j)*tex(i,j+2,1)
      enddo
      enddo
      do i=1,iml
      do j=1,indg+2
         wspsbt(i,j)=(wtbsbt(i,j)+(wtasbt(i,j)-wtbsbt(i,j))*fbcoef)
     &        *tex(i,j+2,1)
      enddo
      enddo
! avoid double counting
        jsn = jndg+5
      else
        jsn = 3
      endif
	  
! north
      if(ip_y .eq. jpe-1) then
      do i=1,iml
      do j=1,indg+2
         htanbt(i,j)=htbnbt(i,j)
     &       +tex(i,jml-1-j,1)*c2dttr*wspnbt(i,j)
      enddo
      enddo
      do i=3,iml-2
      do j=3,jndg+2 !jndg=5
#ifdef BNDTIDE
         htnbt(i,j) = h_tide(i,jml-1-j)
#endif
         hta(i,jml-1-j)=(hta(i,jml-1-j)
     &    +chfht*c2dttr*(htnbt(i,j)-ht(i,jml-1-j)))*tex(i,jml-1-j,1)
      enddo
      enddo
      do i=1,iml
      do j=1,2
#ifdef BNDTIDE
         htanbt(i,j) = h_tide(i,jml-1-j)
#endif
         hta(i,jml-1-j)=htanbt(i,j)*tex(i,jml-1-j,1)
      enddo
      enddo
      do i=1,iml
      do j=1,jndg+2
         wspnbt(i,j)=(wtbnbt(i,j)+(wtanbt(i,j)-wtbnbt(i,j))*fbcoef)
     $        *tex(i,jml-1-j,1)
      enddo
      enddo
! avoid double counting
        jen = jml-jndg-4
      else
        jen = jml-2
      endif

! west
      if(ip_x .eq. 0) then
      do j=1,jml
      do i=1,indg+2
         htawbt(i,j)=htbwbt(i,j)+tex(i+2,j,1)*c2dttr*wspwbt(i,j)
      enddo
      enddo
      do j=jsn,jen
      do i=3,indg+2
#ifdef BNDTIDE
         htwbt(i,j) = h_tide(i+2,j)
#endif
         hta(i+2,j)=(hta(i+2,j)+chfht*c2dttr*(htwbt(i,j)-ht(i+2,j))
     &        )*tex(i+2,j,1)
      enddo
      enddo
      do j=3,jml-2
      do i=1,2
#ifdef BNDTIDE
         htawbt(i,j) = h_tide(i+2,j)
#endif
         hta(i+2,j)=htawbt(i,j)*tex(i+2,j,1)
      enddo
      enddo
      do j=1,jml
      do i=1,indg+2
         wspwbt(i,j)=(wtbwbt(i,j)+(wtawbt(i,j)-wtbwbt(i,j))*fbcoef)
     $        *tex(i+2,j,1)
      enddo
      enddo

      endif
	  
! east
      if(ip_x .eq. ipe-1) then
      do j=1,jml
      do i=1,indg+2
         htaebt(i,j)=htbebt(i,j)
     &       +tex(iml-1-j,j,1)*c2dttr*wspebt(i,j)
      enddo
      enddo
      do j=jsn,jen
      do i=3,indg+2
#ifdef BNDTIDE
         htebt(i,j) = h_tide(iml-1-i,j)
#endif
         hta(iml-1-i,j)=(hta(iml-1-i,j)
     &    +chfht*c2dttr*(htebt(i,j)-ht(iml-1-i,j)))*tex(iml-1-i,j,1)
      enddo
      enddo
      do j=3,jml-2
      do i=1,2
#ifdef BNDTIDE
         htaebt(i,j) = h_tide(iml-1-i,j)
#endif
         hta(iml-1-i,j)=htaebt(i,j)*tex(iml-1-i,j,1)
      enddo
      enddo
      do j=1,jml
      do i=1,indg+2
         wspebt(i,j)=(wtbebt(i,j)+(wtaebt(i,j)-wtbebt(i,j))*fbcoef)
     &        *tex(iml-1-i,j,1)
      enddo
      enddo

      endif
#endif


      call exch_t2d_s1(hta,3)
      call wait_t2d_s1(hta,3)
      call exch_t2d_w1(hta,12)
      call wait_t2d_w1(hta,12)
      call exch_t2d_e1(hta,15) ! debugged by LQH on 09/2010
      call wait_t2d_e1(hta,15) !
      call exch_t2d_n1(hta,16) ! 
      call wait_t2d_n1(hta,16) ! 


      call wait_t2d_e1p(uma,vma,10)
      call wait_t2d_w1p(uma,vma,11)


      if(matsu2.eq.0) then
      do j = 2,jml-1
      do i = 2,iml-1
            umm(i,j) = umm(i,j)+uma(i,j)/dble(ntrop*2+1)
            vmm(i,j) = vmm(i,j)+vma(i,j)/dble(ntrop*2+1)
         enddo
         enddo
      endif

	  
      do j=2,jml-2
      do i=2,iml-2
         hu(i,j)=(ashf(j)*(hta(i,j)+hta(i+1,j))
     &        +anhf(j)*(hta(i,j+1)+hta(i+1,j+1)))
     &        *ex(i,j,1)/areauu(j)
      enddo
      enddo

#ifdef NESTED
! west
      if(ip_x .eq. 0) then
      do j = 1,jml-1
         do i=1,2
          hu(i+2,j)=(ashf(j)*(htawbt(i,j)+htawbt(i+1,j))
     &        +anhf(j)*(htawbt(i,j+1)+htawbt(i+1,j+1)))
     &        *ex(i+2,j,1)/areauu(j)
         enddo
      enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
      do j = 1,jml-1
         do i=1,2
           hu(iml-2-i,j)=(ashf(j)*(htaebt(i,j)+htaebt(i+1,j))
     &        +anhf(j)*(htaebt(i,j+1)+htaebt(i+1,j+1)))
     &        *ex(iml-2-1,j,1)/areauu(j)
         enddo
      enddo
      endif
! south
      if(ip_y .eq. 0) then
      do i=1,iml-1
         do j=1,2
          hu(i,j+2)=(ashf(j+2)*(htasbt(i,j)+htasbt(i+1,j))
     &        +anhf(j+2)*(htasbt(i,j+1)+htasbt(i+1,j+1)))
     &        *ex(i,j+2,1)/areauu(j+2)
         enddo
      enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
      do i=1,iml-1
         do j=1,2
            hu(i,jml-2-j)=(ashf(jml-2-j)*(htanbt(i,j)+htanbt(i+1,j))
     &        +anhf(jml-2-j)*(htanbt(i,j+1)+htanbt(i+1,j+1)))
     &        *ex(i,jml-2-j,1)/areauu(jml-2-j)
         enddo
      enddo
      endif
#endif

      call exch_t2d_n1(hu,4)
      call wait_t2d_n1(hu,4)
      call exch_t2d_e1(hu,13)
      call wait_t2d_e1(hu,13)

#ifdef BNDTIDE
      do i = 1,iml-1
      do j = 1,jml-1
         if(lw(ip_x)+i-1 >= 114 .and. lw(ip_x)+i-1 <= 160 .and.
     &         ls(ip_y)+j-1 >= 335) then
           uma(i,j) = 0.
           vma(i,j) = 0.
           hta(i,j) = 0.
         endif
      enddo
      enddo
#endif

      do j = 2,jml-1
      do i = 2,iml-1
        if(matsu1.ne.1) then
            htb(i,j) = ht(i,j)
            umb(i,j) = um(i,j)
            vmb(i,j) = vm(i,j)
            sfunb(i,j) = sfun(i,j)
            sfvnb(i,j) = sfvn(i,j)
        endif
           ht(i,j) = hta(i,j)
           um(i,j) = uma(i,j)
           vm(i,j) = vma(i,j)
           sfun(i,j)=ex(i,j,1)*uma(i,j)
     &        /(hrr(i,j)+hu(i,j)+(1.d0-ex(i,j,1)))
           sfvn(i,j)=ex(i,j,1)*vma(i,j)
     &        /(hrr(i,j)+hu(i,j)+(1.d0-ex(i,j,1)))
      enddo
      enddo


#ifdef NESTED
! west
      if(ip_x .eq. 0) then
         do j=1,jml
         do i=1,indg+2
         if(matsu1.ne.1) then
            htbwbt(i,j)=htwbt(i,j)
         endif
            htwbt(i,j)=htawbt(i,j)
         enddo
         enddo
      endif
! east
      if(ip_x .eq. ipe-1) then
         do j=1,jml
         do i=1,indg+2
         if(matsu1.ne.1) then
            htbebt(i,j)=htebt(i,j)
         endif
            htebt(i,j)=htaebt(i,j)
         enddo
         enddo
      endif
! south
      if(ip_y .eq. 0) then
         do i=1,iml
         do j=1,jndg+2
         if(matsu1.ne.1) then
            htbsbt(i,j)=htsbt(i,j)
         endif
            htsbt(i,j)=htasbt(i,j)
         enddo
         enddo
      endif
! north
      if(ip_y .eq. jpe-1) then
         do i=1,iml
         do j=1,jndg+2
         if(matsu1.ne.1) then
            htbnbt(i,j)=htnbt(i,j)
         endif
            htnbt(i,j)=htanbt(i,j)
         enddo
         enddo
      endif
#endif


      if(matsu2.eq.1) then
         matsu2 = 0
         go to 1300
      endif

!   ntrop until here
 1200 continue

      do j = 2,jml-1
      do i = 2,iml-1
         um(i,j) = umm(i,j)
         vm(i,j) = vmm(i,j)
      enddo
      enddo

      do j=3,jml-1
      do i=3,iml-1
         wl(i,j,1)=tex(i,j,1)*.25d0*(
     &        dy*((ex(i-1,j,1)+ex(i-1,j-1,1))
     &           *((2.d0-ex(i-1,j-1,1))*um(i-1,j)
     &           +(2.d0-ex(i-1,j,1))*um(i-1,j-1))
     &        -(ex(i,j,1)+ex(i,j-1,1))
     &           *((2.d0-ex(i,j-1,1))*um(i,j)
     &           +(2.d0-ex(i,j,1))*um(i,j-1)))
     &        +dx*((ex(i,j-1,1)+ex(i-1,j-1,1))
     &           *((2.d0-ex(i-1,j-1,1))*vm(i,j-1)
     &           +(2.d0-ex(i,j-1,1))*vm(i-1,j-1))*cs(j-1)
     &        -(ex(i,j,1)+ex(i-1,j,1))
     &           *((2.d0-ex(i-1,j,1))*vm(i,j)
     &           +(2.d0-ex(i,j,1))*vm(i-1,j))*cs(j)) )
      enddo
      enddo

      call exch_t3d_s1(wl,5)
      call wait_t3d_s1(wl,5)
      call exch_t3d_w1(wl,14)

      do j=2,jml-1
      do i=2,iml-1
         sfun(i,j)=ex(i,j,1)*um(i,j)
     &        /(hrr(i,j)+hclu(i,j)+(1.d0-ex(i,j,1)))
         sfvn(i,j)=ex(i,j,1)*vm(i,j)
     &        /(hrr(i,j)+hclu(i,j)+(1.d0-ex(i,j,1)))
      enddo
      enddo

  
!   form u/v
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
         u(i,j,k)=ex(i,j,k)*(u(i,j,k)+sfun(i,j))
         v(i,j,k)=ex(i,j,k)*(v(i,j,k)+sfvn(i,j))
#ifdef DEBUG1
#ifdef GL11M
         if(abs(u(i,j,k)).ge.300.) then
            write(*,*) 'u# (trop): ',nkai,ip,i,j,
     &  lw(ip_x)+i-1,ls(ip_y)+j-1,k,u(i,j,k)
         endif
         if(abs(v(i,j,k)).ge.300.) then
            write(*,*) 'v# (trop): ',nkai,ip,i,j,
     &  lw(ip_x)+i-1,ls(ip_y)+j-1,k,v(i,j,k)
         endif
#endif
         uvlimit=350.d0
         if(abs(u(i,j,k)).ge.uvlimit) then
            write(*,*) 'u# (trop): ',
     &         nkai,lw(ip_x)+i-1,ls(ip_y)+j-1,k,u(i,j,k),sfun(i,j),
     &           sfu(i,j)
         endif
         if(abs(v(i,j,k)).ge.uvlimit) then
            write(*,*) 'v# (trop): ',
     &         nkai,lw(ip_x)+i-1,ls(ip_y)+j-1,k,v(i,j,k),sfvn(i,j),
     &           sfv(i,j)
         endif
#endif

      enddo
      enddo
      enddo

      call wait_t3d_w1(wl,14)

      return
      end