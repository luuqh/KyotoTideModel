ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     isopycnal diffusion with eddy parameterization                c
c     by Gent and McWillams (1990)                                  c
c     output horizontal and vertical diffusion for ta,sa            c
c     request for tracli                                            c
c                                                                   c
c     $Id: gm_dif.F90 2 2007-11-15 13:07:08Z ishikawa $
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gm_dif

      use param
      use mod_mpi
      implicit none

#include "common.h"

      integer i,j,k
      real*8 druws,tuuw,tvus,suuw,svus,sxue,sxuw,syun,syus,
     &     drues,tuue,suue,druwn,tvun,svun,
     &     drdws,tudw,tvds,sudw,svds,sxde,sxdw,sydn,syds,
     &     drdes,tude,sude,drdwn,tvdn,svdn,
     &     vdf,tvt2,svt2,sxw,sxe,sys,syn,azz,drho,
     &     fchdf,tvt2mx,svt2mx,azzt,azzs,
     &     dtdx1(im,jm,km+1),dsdx1(im,jm,km+1),
     &     dtdy1(im,jm,km+1),dsdy1(im,jm,km+1),
     &     dtdz1(im,jm,km+1),dsdz1(im,jm,km+1)
c
c
ccc   facter to lateral mixing for calm spin-up
      fchdf=0.d0
#ifdef PC68M
      fchdf=dmax1( dble(nxmonth-nkai)/dble(nxmonth),0.d0 )
#endif
c
ccc   isopycnal diffusion with gm-transport
c
      do k = 1,km+1
      do j = 2,jml-1
      do i = 2,iml-1
      dtdx1(i,j,k) = tex(i+1,j,k)*tex(i,j,k)
     &       *(tb(i+1,j,k)-tb(i,j,k))*dxr*cstr(j)
      dsdx1(i,j,k) = tex(i+1,j,k)*tex(i,j,k)
     &       *(sb(i+1,j,k)-sb(i,j,k))*dxr*cstr(j)
      dtdy1(i,j,k) = tex(i,j+1,k)*tex(i,j,k)
     &       *(tb(i,j+1,k)-tb(i,j,k))*dyr
      dsdy1(i,j,k) = tex(i,j+1,k)*tex(i,j,k)
     &       *(sb(i,j+1,k)-sb(i,j,k))*dyr
      dtdz1(i,j,k) = tex(i,j,k)*(tb(i,j,k-1)-tb(i,j,k))
     &        *2./(dzt(i,j,k-1)+dzt(i,j,k)+1.-tex(i,j,k-1))
      dsdz1(i,j,k) = tex(i,j,k)*(sb(i,j,k-1)-sb(i,j,k))
     &        *2./(dzt(i,j,k-1)+dzt(i,j,k)+1.-tex(i,j,k-1))
      enddo
      enddo
      enddo

      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,k).eq.1.)then
c
      druws = -tex(i,j,k)*(rhodf(i,j,k-1)-rho(i,j,k))
     &        *2.d0/(dzt(i,j,k-1)+dzt(i,j,k)
     &        +(1.d0-tex(i,j,k-1))*(1.d0-tex(i,j,k)) )
c
      if(druws.le.drmax) then
         tuuw = hdtsmin*(dzu(i,j-1,k)+dzu(i,j,k))*dtdx1(i,j,k)
         tvus = hdtsmin*(dzu(i-1,j,k)+dzu(i,j,k))*dtdy1(i,j,k)
         suuw = hdtsmin*(dzu(i,j-1,k)+dzu(i,j,k))*dsdx1(i,j,k)
         svus = hdtsmin*(dzu(i-1,j,k)+dzu(i,j,k))*dsdy1(i,j,k)
      else
         sxue = tex(i+1,j,k)*tex(i,j,k)*
     &        (rho(i+1,j,k)-rho(i,j,k))*dxr*cstr(j)/druws
         sxuw = tex(i,j,k)*tex(i-1,j,k)*
     &        (rho(i,j,k)-rho(i-1,j,k))*dxr*cstr(j)/druws
         syun = tex(i,j+1,k)*tex(i,j,k)*
     &        (rho(i,j+1,k)-rho(i,j,k))*dyr/druws
         syus = tex(i,j,k)*tex(i,j-1,k)*
     &        (rho(i,j,k)-rho(i,j-1,k))*dyr/druws
         sxue=dsign(dmin1(dabs(sxue),sxymax),sxue)
         sxuw=dsign(dmin1(dabs(sxuw),sxymax),sxuw)
         syun=dsign(dmin1(dabs(syun),sxymax),syun)
         syus=dsign(dmin1(dabs(syus),sxymax),syus)
#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*)  'druws',druws,sxue,sxuw,syun,syus
      endif
#endif
         tuuw = aip*(dzu(i,j-1,k)/(1.+sxue*sxue+syus*syus)*(
     &        (1.+syus*syus+ar*sxue*sxue)*dtdx1(i,j,k)
     &        +(ar-1.)*sxue*syus         *dtdy1(i,j-1,k)
     &        +(1.-ar)*sxue              *dtdz1(i,j,k) )
     &        +dzu(i,j,k)/(1.+sxue*sxue+syun*syun)*(
     &        (1.+syun*syun+ar*sxue*sxue)*dtdx1(i,j,k)
     &        +(ar-1.)*sxue*syun         *dtdy1(i,j,k)
     &        +(1.-ar)*sxue              *dtdz1(i,j,k) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxue
     &                                   *dtdz1(i,j,k)
c
         suuw = aip*(dzu(i,j-1,k)/(1.+sxue*sxue+syus*syus)*(
     &        (1.+syus*syus+ar*sxue*sxue)*dsdx1(i,j,k)
     &        +(ar-1.)*sxue*syus*dsdy1(i,j-1,k)
     &        +(1.-ar)*sxue*dsdz1(i,j,k) )
     &        +dzu(i,j,k)/(1.+sxue*sxue+syun*syun)*(
     &        (1.+syun*syun+ar*sxue*sxue)*dsdx1(i,j,k)
     &        +(ar-1.)*sxue*syun*dsdy1(i,j,k)
     &        +(1.-ar)*sxue*dsdz1(i,j,k) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxue
     &        *dsdz1(i,j,k)

         tvus = aip*(dzu(i-1,j,k)/(1.+sxuw*sxuw+syun*syun)*(
     &        +(ar-1.)*sxuw*syun*dtdx1(i-1,j,k)
     &        +(1.+sxuw*sxuw+ar*syun*syun)*dtdy1(i,j,k)
     &        +(1.-ar)*syun*dtdz1(i,j,k) )
     &        +dzu(i,j,k)/(1.+sxue*sxue+syun*syun)*(
     &        +(ar-1.)*sxue*syun*dtdx1(i,j,k)
     &        +(1.+sxue*sxue+ar*syun*syun)*dtdy1(i,j,k)
     &        +(1.-ar)*syun*dtdz1(i,j,k) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*syun
     &        *dtdz1(i,j,k)
c
         svus = aip*(dzu(i-1,j,k)/(1.+sxuw*sxuw+syun*syun)*(
     &        (ar-1.)*sxuw*syun*dsdx1(i-1,j,k)
     &        +(1.+sxuw*sxuw+ar*syun*syun)*dsdy1(i,j,k)
     &        +(1.-ar)*syun*dsdz1(i,j,k) )
     &        +dzu(i,j,k)/(1.+sxue*sxue+syun*syun)*(
     &        (ar-1.)*sxue*syun*dsdx1(i,j,k)
     &        +(1.+sxue*sxue+ar*syun*syun)*dsdy1(i,j,k)
     &        +(1.-ar)*syun*dsdz1(i,j,k) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*syun
     &        *dsdz1(i,j,k)
      endif
c
      drues = -tex(i+1,j,k)*(rhodf(i+1,j,k-1)-rho(i+1,j,k))
     &     *2./(dzt(i+1,j,k-1)+dzt(i+1,j,k)
     &     +(1.-tex(i+1,j,k-1))*(1.-tex(i+1,j,k)) )
c
      if(drues.le.drmax) then
         tuue = hdtsmin*(dzu(i,j-1,k)+dzu(i,j,k))*dtdx1(i,j,k)
         suue = hdtsmin*(dzu(i,j-1,k)+dzu(i,j,k))*dsdx1(i,j,k)
      else
         sxuw = tex(i+1,j,k)*tex(i,j,k)*
     &        (rho(i+1,j,k)-rho(i,j,k))*dxr*cstr(j)/drues
         syun = tex(i+1,j+1,k)*tex(i+1,j,k)*
     &        (rho(i+1,j+1,k)-rho(i+1,j,k))*dyr/drues
         syus = tex(i+1,j,k)*tex(i+1,j-1,k)*
     &        (rho(i+1,j,k)-rho(i+1,j-1,k))*dyr/drues
         sxuw=dsign(dmin1(dabs(sxuw),sxymax),sxuw)
         syun=dsign(dmin1(dabs(syun),sxymax),syun)
         syus=dsign(dmin1(dabs(syus),sxymax),syus)
#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*)  'drues',drues,sxuw,syun,syus
           write(*,*) 'syun',rho(i+1,j+1,k),rho(i+1,j,k),
     &       t(i+1,j+1,k),t(i+1,j,k),s(i+1,j+1,k),s(i+1,j,k)
      endif
#endif
         tuue = aip*(dzu(i,j-1,k)/(1.+sxuw*sxuw+syus*syus)*(
     &        (1.+syus*syus+ar*sxuw*sxuw)*dtdx1(i,j,k)
     &        +(ar-1.)*sxuw*syus*dtdy1(i+1,j-1,k)
     &        +(1.-ar)*sxuw*dtdz1(i+1,j,k) )
     &        +dzu(i,j,k)/(1.+sxuw*sxuw+syun*syun)*(
     &        (1.+syun*syun+ar*sxuw*sxuw)*dtdx1(i,j,k)
     &        +(ar-1.)*sxuw*syun*dtdy1(i+1,j,k)
     &        +(1.-ar)*sxuw*dtdz1(i+1,j,k) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxuw
     &        *dtdz1(i+1,j,k)
c
         suue = aip*(dzu(i,j-1,k)/(1.+sxuw*sxuw+syus*syus)*(
     &        (1.+syus*syus+ar*sxuw*sxuw)*dsdx1(i,j,k)
     &        +(ar-1.)*sxuw*syus*dsdy1(i+1,j-1,k)
     &        +(1.-ar)*sxuw*dsdz1(i+1,j,k) )
     &        +dzu(i,j,k)/(1.+sxuw*sxuw+syun*syun)*(
     &        (1.+syun*syun+ar*sxuw*sxuw)*dsdx1(i,j,k)
     &        +(ar-1.)*sxuw*syun*dsdy1(i+1,j,k)
     &        +(1.-ar)*sxuw*dsdz1(i+1,j,k) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxuw
     &        *dsdz1(i+1,j,k)
      endif
c
      druwn = -tex(i,j+1,k)*(rhodf(i,j+1,k-1)-rho(i,j+1,k))
     &     *2./(dzt(i,j+1,k-1)+dzt(i,j+1,k)
     &     +(1.-tex(i,j+1,k-1))*(1.-tex(i,j+1,k)) )
c
      if(druwn.le.drmax) then
         tvun = hdtsmin*(dzu(i-1,j,k)+dzu(i,j,k))*dtdy1(i,j,k)
         svun = hdtsmin*(dzu(i-1,j,k)+dzu(i,j,k))*dsdy1(i,j,k)
      else
         sxue = tex(i+1,j+1,k)*tex(i,j+1,k)*
     &        (rho(i+1,j+1,k)-rho(i,j+1,k))*dxr*cstr(j+1)/druwn
         sxuw = tex(i,j+1,k)*tex(i-1,j+1,k)*
     &        (rho(i,j+1,k)-rho(i-1,j+1,k))*dxr*cstr(j+1)/druwn
         syus = tex(i,j+1,k)*tex(i,j,k)*
     &        (rho(i,j+1,k)-rho(i,j,k))*dyr/druwn
         sxue=dsign(dmin1(dabs(sxue),sxymax),sxue)
         sxuw=dsign(dmin1(dabs(sxuw),sxymax),sxuw)
         syus=dsign(dmin1(dabs(syus),sxymax),syus)
#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*)  'druwn',druwn,sxue,sxuw,syus
           write(*,*) 'sxue',rho(i+1,j+1,k),rho(i,j+1,k),
     &   t(i+1,j+1,k),t(i,j+1,k),s(i+1,j+1,k),s(i,j+1,k)
      endif
#endif

         tvun = aip*(dzu(i-1,j,k)/(1.+sxuw*sxuw+syus*syus)*(
     &        +(ar-1.)*sxuw*syus*dtdx1(i-1,j+1,k)
     &        +(1.+sxuw*sxuw+ar*syus*syus)*dtdy1(i,j,k)
     &        +(1.-ar)*syus*dtdz1(i,j+1,k) )
     &        +dzu(i,j,k)/(1.+sxue*sxue+syus*syus)*(
     &        +(ar-1.)*sxue*syus*dtdx1(i,j+1,k)
     &        +(1.+sxue*sxue+ar*syus*syus)*dtdy1(i,j,k)
     &        +(1.-ar)*syus*dtdz1(i,j+1,k) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*syus
     &        *dtdz1(i,j+1,k)
c
         svun = aip*(dzu(i-1,j,k)/(1.+sxuw*sxuw+syus*syus)*(
     &        +(ar-1.)*sxuw*syus*dsdx1(i-1,j+1,k)
     &        +(1.+sxuw*sxuw+ar*syus*syus)*dsdy1(i,j,k)
     &        +(1.-ar)*syus*dsdz1(i,j+1,k))
     &        +dzu(i,j,k)/(1.+sxue*sxue+syus*syus)*(
     &        +(ar-1.)*sxue*syus*dsdx1(i,j+1,k)
     &        +(1.+sxue*sxue+ar*syus*syus)*dsdy1(i,j,k)
     &        +(1.-ar)*syus*dsdz1(i,j+1,k) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*syus
     &        *dsdz1(i,j+1,k)
      endif
c
      drdws = -tex(i,j,k)*(rho(i,j,k)-rhouf(i,j,k+1))
     &     *2./(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
c
      if(drdws.le.drmax) then
         tudw = (tex(i,j,k+1)*hdtsmin+(1.-tex(i,j,k+1))*aip)
     &        *(dzu(i,j-1,k)+dzu(i,j,k))*dtdx1(i,j,k)
         tvds = (tex(i,j,k+1)*hdtsmin+(1.-tex(i,j,k+1))*aip)
     &        *(dzu(i-1,j,k)+dzu(i,j,k))*dtdy1(i,j,k)
         sudw = (tex(i,j,k+1)*hdtsmin+(1.-tex(i,j,k+1))*aip)
     &        *(dzu(i,j-1,k)+dzu(i,j,k))*dsdx1(i,j,k)
         svds = (tex(i,j,k+1)*hdtsmin+(1.-tex(i,j,k+1))*aip)
     &        *(dzu(i-1,j,k)+dzu(i,j,k))*dsdy1(i,j,k)
      else
         sxde = tex(i+1,j,k)*tex(i,j,k+1)*
     &        (rho(i+1,j,k)-rho(i,j,k))*dxr*cstr(j)/drdws
         sxdw = tex(i,j,k+1)*tex(i-1,j,k)*
     &        (rho(i,j,k)-rho(i-1,j,k))*dxr*cstr(j)/drdws
         sydn = tex(i,j+1,k)*tex(i,j,k+1)*
     &        (rho(i,j+1,k)-rho(i,j,k))*dyr/drdws
         syds = tex(i,j,k+1)*tex(i,j-1,k)*
     &        (rho(i,j,k)-rho(i,j-1,k))*dyr/drdws
         sxde=dsign(dmin1(dabs(sxde),sxymax),sxde)
         sxdw=dsign(dmin1(dabs(sxdw),sxymax),sxdw)
         sydn=dsign(dmin1(dabs(sydn),sxymax),sydn)
         syds=dsign(dmin1(dabs(syds),sxymax),syds)
#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*)  'drdws',drdws,sxde,sxdw,sydn,syds
      endif
#endif

         tudw = aip*(dzu(i,j-1,k)/(1.+sxde*sxde+syds*syds)*(
     &        (1.+syds*syds+ar*sxde*sxde)*dtdx1(i,j,k)
     &        +(ar-1.)*sxde*syds*dtdy1(i,j-1,k)
     &        +(1.-ar)*sxde*dtdz1(i,j,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+sydn*sydn)*(
     &        (1.+sydn*sydn+ar*sxde*sxde)*dtdx1(i,j,k)
     &        +(ar-1.)*sxde*sydn*dtdy1(i,j,k)
     &        +(1.-ar)*sxde*dtdz1(i,j,k+1) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxde
     &        *dtdz1(i,j,k+1)
c
         sudw = aip*(dzu(i,j-1,k)/(1.+sxde*sxde+syds*syds)*(
     &        (1.+syds*syds+ar*sxde*sxde)*dsdx1(i,j,k)
     &        +(ar-1.)*sxde*syds*dsdy1(i,j-1,k)
     &        +(1.-ar)*sxde*dsdz1(i,j,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+sydn*sydn)*(
     &        (1.+sydn*sydn+ar*sxde*sxde)*dsdx1(i,j,k)
     &        +(ar-1.)*sxde*sydn*dsdy1(i,j,k)
     &        +(1.-ar)*sxde*dsdz1(i,j,k+1) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxde
     &        *dsdz1(i,j,k+1)
c
         tvds = aip*(dzu(i-1,j,k)/(1.+sxdw*sxdw+sydn*sydn)*(
     &        +(ar-1.)*sxdw*sydn*dtdx1(i-1,j,k)
     &        +(1.+sxdw*sxdw+ar*sydn*sydn)*dtdy1(i,j,k)
     &        +(1.-ar)*sydn*dtdz1(i,j,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+sydn*sydn)*(
     &        +(ar-1.)*sxde*sydn*dtdx1(i,j,k)
     &        +(1.+sxde*sxde+ar*sydn*sydn)*dtdy1(i,j,k)
     &        +(1.-ar)*sydn*dtdz1(i,j,k+1) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*sydn
     &        *dtdz1(i,j,k+1)
c
         svds = aip*(dzu(i-1,j,k)/(1.+sxdw*sxdw+sydn*sydn)*(
     &        +(ar-1.)*sxdw*sydn*dsdx1(i-1,j,k)
     &        +(1.+sxdw*sxdw+ar*sydn*sydn)*dsdy1(i,j,k)
     &        +(1.-ar)*sydn*dsdz1(i,j,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+sydn*sydn)*(
     &        +(ar-1.)*sxde*sydn*dsdx1(i,j,k)
     &        +(1.+sxde*sxde+ar*sydn*sydn)*dsdy1(i,j,k)
     &        +(1.-ar)*sydn*dsdz1(i,j,k+1) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*sydn
     &        *dsdz1(i,j,k+1)
      endif
c
      drdes = -tex(i+1,j,k)*(rho(i+1,j,k)-rhouf(i+1,j,k+1))
     &     *2./(dzt(i+1,j,k)+dzt(i+1,j,k+1)+1.-tex(i+1,j,k))
c
      if(drdes.le.drmax) then
         tude = (tex(i+1,j,k+1)*hdtsmin+(1.-tex(i+1,j,k+1))*aip)
     &        *(dzu(i,j-1,k)+dzu(i,j,k))*dtdx1(i,j,k)
         sude = (tex(i+1,j,k+1)*hdtsmin+(1.-tex(i+1,j,k+1))*aip)
     &        *(dzu(i,j-1,k)+dzu(i,j,k))*dsdx1(i,j,k)
      else
         sxde = tex(i,j,k)*tex(i+1,j,k+1)*
     &        (rho(i+1,j,k)-rho(i,j,k))*dxr*cstr(j)/drdes
         sydn = tex(i+1,j+1,k)*tex(i+1,j,k+1)*
     &        (rho(i+1,j+1,k)-rho(i+1,j,k))*dyr/drdes
         syds = tex(i+1,j-1,k)*tex(i+1,j,k+1)*
     &        (rho(i+1,j,k)-rho(i+1,j-1,k))*dyr/drdes
         sxde=dsign(dmin1(dabs(sxde),sxymax),sxde)
         sydn=dsign(dmin1(dabs(sydn),sxymax),sydn)
         syds=dsign(dmin1(dabs(syds),sxymax),syds)
#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*)  'drdes',drdes,sxde,sydn,syds
      endif
#endif
         tude = aip*(dzu(i,j-1,k)/(1.+sxde*sxde+syds*syds)*(
     &        (1.+syds*syds+ar*sxde*sxde)*dtdx1(i,j,k)
     &        +(ar-1.)*sxde*syds*dtdy1(i+1,j-1,k)
     &        +(1.-ar)*sxde*dtdz1(i+1,j,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+sydn*sydn)*(
     &        (1.+sydn*sydn+ar*sxde*sxde)*dtdx1(i,j,k)
     &        +(ar-1.)*sxde*sydn*dtdy1(i+1,j,k)
     &        +(1.-ar)*sxde*dtdz1(i+1,j,k+1) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxde
     &        *dtdz1(i+1,j,k+1)
c
         sude = aip*(dzu(i,j-1,k)/(1.+sxde*sxde+syds*syds)*(
     &        (1.+syds*syds+ar*sxde*sxde)*dsdx1(i,j,k)
     &        +(ar-1.)*sxde*syds*dsdy1(i+1,j-1,k)
     &        +(1.-ar)*sxde*dsdz1(i+1,j,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+sydn*sydn)*(
     &        (1.+sydn*sydn+ar*sxde*sxde)*dsdx1(i,j,k)
     &        +(ar-1.)*sxde*sydn*dsdy1(i+1,j,k)
     &        +(1.-ar)*sxde*dsdz1(i+1,j,k+1) ))
     &        -agm*(dzu(i,j-1,k)+dzu(i,j,k))*sxde
     &        *dsdz1(i+1,j,k+1)
      endif
c
      drdwn = -tex(i,j+1,k)*(rho(i,j+1,k)-rhouf(i,j+1,k+1))
     &     *2./(dzt(i,j+1,k)+dzt(i,j+1,k+1)+1.-tex(i,j+1,k))
c
      if(drdwn.le.drmax) then
         tvdn = (tex(i,j+1,k+1)*hdtsmin+(1.-tex(i,j+1,k+1))*aip)
     &        *(dzu(i-1,j,k)+dzu(i,j,k))*dtdy1(i,j,k)
         svdn = (tex(i,j+1,k+1)*hdtsmin+(1.-tex(i,j+1,k+1))*aip)
     &        *(dzu(i-1,j,k)+dzu(i,j,k))*dsdy1(i,j,k)
      else
         sxde = tex(i+1,j+1,k)*tex(i,j+1,k+1)*
     &        (rho(i+1,j+1,k)-rho(i,j+1,k))*dxr*cstr(j+1)/drdwn
         sxdw = tex(i,j+1,k+1)*tex(i-1,j+1,k)*
     &        (rho(i,j+1,k)-rho(i-1,j+1,k))*dxr*cstr(j+1)/drdwn
         syds = tex(i,j+1,k+1)*tex(i,j,k)*
     &        (rho(i,j+1,k)-rho(i,j,k))*dyr/drdwn
         sxde=dsign(dmin1(dabs(sxde),sxymax),sxde)
         sxdw=dsign(dmin1(dabs(sxdw),sxymax),sxdw)
         syds=dsign(dmin1(dabs(syds),sxymax),syds)

#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*)  'drdwn',sxde,sxdw,syds,drdwn
      endif
#endif

         tvdn = aip*(dzu(i-1,j,k)/(1.+sxdw*sxdw+syds*syds)*(
     &        (ar-1.)*sxdw*syds*dtdx1(i-1,j+1,k)
     &        +(1.+sxdw*sxdw+ar*syds*syds)*dtdy1(i,j,k)
     &        +(1.-ar)*syds*dtdz1(i,j+1,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+syds*syds)*(
     &        (ar-1.)*sxde*syds*dtdx1(i,j+1,k)
     &        +(1.+sxde*sxde+ar*syds*syds)*dtdy1(i,j,k)
     &        +(1.-ar)*syds*dtdz1(i,j+1,k+1) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*syds
     &        *dtdz1(i,j+1,k+1)
c
         svdn = aip*(dzu(i-1,j,k)/(1.+sxdw*sxdw+syds*syds)*(
     &        (ar-1.)*sxdw*syds*dsdx1(i-1,j+1,k)
     &        +(1.+sxdw*sxdw+ar*syds*syds)*dsdy1(i,j,k)
     &        +(1.-ar)*syds*dsdz1(i,j+1,k+1) )
     &        +dzu(i,j,k)/(1.+sxde*sxde+syds*syds)*(
     &        (ar-1.)*sxde*syds*dsdx1(i,j+1,k)
     &        +(1.+sxde*sxde+ar*syds*syds)*dsdy1(i,j,k)
     &        +(1.-ar)*syds*dsdz1(i,j+1,k+1) ))
     &        -agm*(dzu(i-1,j,k)+dzu(i,j,k))*syds
     &        *dsdz1(i,j+1,k+1)
      endif
c
      tudf(i,j,k) = -0.125*dy*(tuuw+tuue+tudw+tude)
     &     *tex(i,j,k)*tex(i+1,j,k)
      tvdf(i,j,k) = -0.125*dx*cs(j)*(tvus+tvun+tvds+tvdn)
     &     *tex(i,j,k)*tex(i,j+1,k)
      sudf(i,j,k) = -0.125*dy*(suuw+suue+sudw+sude)
     &     *tex(i,j,k)*tex(i+1,j,k)
      svdf(i,j,k) = -0.125*dx*cs(j)*(svus+svun+svds+svdn)
     &     *tex(i,j,k)*tex(i,j+1,k)
#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*) nkai,'tudf',tudf(i,j,k),tuuw,tuue,tudw,tude
           write(*,*) 'tvdf',tvdf(i,j,k),tvus,tvun,tvds,tvdn
           write(*,*) 'sudf',sudf(i,j,k),suuw,suue,sudw,sude
           write(*,*) 'svdf',svdf(i,j,k),svus,svun,svds,svdn
           write(*,*)  
      endif
#endif
      else
c
      tudf(i,j,k)=0.d0
      tvdf(i,j,k)=0.d0
      sudf(i,j,k)=0.d0
      svdf(i,j,k)=0.d0
c
      endif
c
ccc   lateral diffusion at bottom
      if(tex(i,j,k).eq.1..and.tex(i,j,k+1).eq.0.)then

      tudf(i,j,k)=tudf(i,j,k)
     $        -0.5d0*dy*hdts_btm*(dzu(i,j,k)+dzu(i,j-1,k))
     $        *dtdx1(i,j,k)
      sudf(i,j,k)=sudf(i,j,k)
     $     -0.5d0*dy*hdts_btm*(dzu(i,j,k)+dzu(i,j-1,k))
     $     *dsdx1(i,j,k)
      tvdf(i,j,k)=tvdf(i,j,k)
     $     -.5d0*dx*hdts_btm*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
     $     *dtdy1(i,j,k)
      svdf(i,j,k)=svdf(i,j,k)
     $     -.5d0*dx*hdts_btm*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
     $     *dsdy1(i,j,k)
      endif
c
ccc   lateral diffusion at 1st layer
      if(k.eq.1.and.tex(i,j,k).eq.1.)then
      tudf(i,j,k)=tudf(i,j,k)
     $        -0.5d0*dy*hdts_sfc*(dzu(i,j,k)+dzu(i,j-1,k))
     $        *dtdx1(i,j,k)
      sudf(i,j,k)=sudf(i,j,k)
     $     -0.5d0*dy*hdts_sfc*(dzu(i,j,k)+dzu(i,j-1,k))
     $     *dsdx1(i,j,k)
      tvdf(i,j,k)=tvdf(i,j,k)
     $     -.5d0*dx*hdts_sfc*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
     $     *dtdy1(i,j,k)
      svdf(i,j,k)=svdf(i,j,k)
     $     -.5d0*dx*hdts_sfc*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
     $     *dsdy1(i,j,k)
      endif

ccc   lateral diffusion at northern most
!      if(ls(ip_y)+j-1>561 .and. tex(i,j,k).eq.1.)then
!      tudf(i,j,k)=tudf(i,j,k)
!     $        -0.5d0*dy*hdts_btm*(dzu(i,j,k)+dzu(i,j-1,k))
!     $        *dtdx1(i,j,k)
!      sudf(i,j,k)=sudf(i,j,k)
!     $     -0.5d0*dy*hdts_btm*(dzu(i,j,k)+dzu(i,j-1,k))
!     $     *dsdx1(i,j,k)
!      tvdf(i,j,k)=tvdf(i,j,k)
!     $     -.5d0*dx*hdts_btm*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
!     $     *dtdy1(i,j,k)
!      svdf(i,j,k)=svdf(i,j,k)
!     $     -.5d0*dx*hdts_btm*(dzu(i-1,j,k)+dzu(i,j,k))*cs(j)
!     $     *dsdy1(i,j,k)
!      endif

#ifdef DEBUG1
      if(nkai>482085 .and. nkai < 482100 .and. ip==59 .and.
     &   i==32 .and. j==50 .and. k==22) then
           write(*,*) nkai,'tudf_last',tudf(i,j,k),tvdf(i,j,k),
     &     sudf(i,j,k),svdf(i,j,k)
      endif
#endif

      enddo
      enddo
      enddo

      call exch_t3d_s1p(tvdf,svdf,1,2)
      call exch_t3d_w1p(tudf,sudf,3,4)

      do k = 1,km-1
      do j = 3,jml-2
      do i = 3,iml-2
      if(tex(i,j,k+1).eq.1.)then
c
      tvt2=0.
      svt2=0.
      sxw=0.
      sxe=0.
      sys=0.
      syn=0.
c
      if(vdts(i,j,k+1).ge.vdtsmx) then
         vdts(i,j,k+1)=vdtsmx*tex(i,j,k+1)
         vddt(i,j,k+1)=vdtsmx*tex(i,j,k+1)
         vdds(i,j,k+1)=vdtsmx*tex(i,j,k+1)
         vdf=vdts(i,j,k+1)*areat(i,j,k+1)*2.d0
     $        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
         tvt2=-(tb(i,j,k)-tb(i,j,k+1))*vdf
         svt2=-(sb(i,j,k)-sb(i,j,k+1))*vdf
      else
         drho = -tex(i,j,k+1)*(rhodh(i,j,k)-rhouh(i,j,k+1))
     &        *2./(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
         if(drho.le.drmax) then
            vdts(i,j,k+1)=vdtsmx*tex(i,j,k+1)
            vddt(i,j,k+1)=vdtsmx*tex(i,j,k+1)
            vdds(i,j,k+1)=vdtsmx*tex(i,j,k+1)
            vdf=vdts(i,j,k+1)*areat(i,j,k+1)*2.d0
     $           /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
            tvt2=-(tb(i,j,k)-tb(i,j,k+1))*vdf
            svt2=-(sb(i,j,k)-sb(i,j,k+1))*vdf
         else
            sxw=(ex(i-1,j-1,k+1)+ex(i-1,j,k+1)
     &           -ex(i-1,j-1,k+1)*ex(i-1,j,k+1))
     &           *(ex(i,j-1,k+1)+ex(i,j,k+1)
     &           -ex(i,j-1,k+1)*ex(i,j,k+1))
     &           *((rhodh(i,j,k)+rhouh(i,j,k+1))
     &           -(rhodh(i-1,j,k)+rhouh(i-1,j,k+1)))
     &           *dx2r*cstr(j)/drho
            sxe=(ex(i,j-1,k+1)+ex(i,j,k+1)
     &           -ex(i,j-1,k+1)*ex(i,j,k+1))
     &           *(ex(i-1,j-1,k+1)+ex(i-1,j,k+1)
     &           -ex(i-1,j-1,k+1)*ex(i-1,j,k+1))
     &           *((rhodh(i+1,j,k)+rhouh(i+1,j,k+1))
     &           -(rhodh(i,j,k)+rhouh(i,j,k+1)))*dx2r*cstr(j)/drho
            sys=(ex(i-1,j-1,k+1)+ex(i,j-1,k+1)
     &           -ex(i-1,j-1,k+1)*ex(i,j-1,k+1))
     &           *(ex(i-1,j,k+1)+ex(i,j,k+1)
     &           -ex(i-1,j,k+1)*ex(i,j,k+1))
     &           *((rhodh(i,j,k)+rhouh(i,j,k+1))
     &           -(rhodh(i,j-1,k)+rhouh(i,j-1,k+1)))*dy2r/drho
            syn=(ex(i-1,j,k+1)+ex(i,j,k+1)
     &           -ex(i-1,j,k+1)*ex(i,j,k+1))
     &           *(ex(i-1,j-1,k+1)+ex(i,j-1,k+1)
     &           -ex(i-1,j-1,k+1)*ex(i,j-1,k+1))
     &           *((rhodh(i,j+1,k)+rhouh(i,j+1,k+1))
     &           -(rhodh(i,j,k)+rhouh(i,j,k+1)))*dy2r/drho
            sxw=dsign(dmin1(dabs(sxw),sxymax),sxw)
            sxe=dsign(dmin1(dabs(sxe),sxymax),sxe)
            sys=dsign(dmin1(dabs(sys),sxymax),sys)
            syn=dsign(dmin1(dabs(syn),sxymax),syn)
c
            azz = vdts(i,j,k+1)-vdtsmn + aip*0.25*(
     &            (ar+sxw*sxw+sys*sys)/(1.+sxw*sxw+sys*sys)
     &           +(ar+sxw*sxw+syn*syn)/(1.+sxw*sxw+syn*syn)
     &           +(ar+sxe*sxe+sys*sys)/(1.+sxe*sxe+sys*sys)
     &           +(ar+sxe*sxe+syn*syn)/(1.+sxe*sxe+syn*syn) )
            azzt = vddt(i,j,k+1)-vdtsmn + aip*0.25*(
     &            (ar+sxw*sxw+sys*sys)/(1.+sxw*sxw+sys*sys)
     &           +(ar+sxw*sxw+syn*syn)/(1.+sxw*sxw+syn*syn)
     &           +(ar+sxe*sxe+sys*sys)/(1.+sxe*sxe+sys*sys)
     &           +(ar+sxe*sxe+syn*syn)/(1.+sxe*sxe+syn*syn) )
            azzs = vdds(i,j,k+1)-vdtsmn + aip*0.25*(
     &            (ar+sxw*sxw+sys*sys)/(1.+sxw*sxw+sys*sys)
     &           +(ar+sxw*sxw+syn*syn)/(1.+sxw*sxw+syn*syn)
     &           +(ar+sxe*sxe+sys*sys)/(1.+sxe*sxe+sys*sys)
     &           +(ar+sxe*sxe+syn*syn)/(1.+sxe*sxe+syn*syn) )
            if(azz.ge.vdtsmx) then
               vdts(i,j,k+1)=vdtsmx*tex(i,j,k+1)
               vddt(i,j,k+1)=vdtsmx*tex(i,j,k+1)
               vdds(i,j,k+1)=vdtsmx*tex(i,j,k+1)
               vdf=vdts(i,j,k+1)*areat(i,j,k+1)*2.d0
     $              /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
               tvt2=-(tb(i,j,k)-tb(i,j,k+1))*vdf
               svt2=-(sb(i,j,k)-sb(i,j,k+1))*vdf
            else
c
c         vdts(i,j,k+1)=dmax1(vdts(i,j,k+1),azz,vdtsmn)*tex(i,j,k+1)
c         vddt(i,j,k+1)=dmax1(vddt(i,j,k+1),azz,vdtsmn)*tex(i,j,k+1)
c         vdds(i,j,k+1)=dmax1(vdds(i,j,k+1),azz,vdtsmn)*tex(i,j,k+1)
         vdts(i,j,k+1)=dmax1(azz,vdtsmn)
         vddt(i,j,k+1)=dmax1(azzt,vdtsmn)
         vdds(i,j,k+1)=dmax1(azzs,vdtsmn)
         vdts(i,j,k+1)=dmin1(vdts(i,j,k+1),vdtsmx)*tex(i,j,k+1)
         vddt(i,j,k+1)=dmin1(vddt(i,j,k+1),vdtsmx)*tex(i,j,k+1)
         vdds(i,j,k+1)=dmin1(vdds(i,j,k+1),vdtsmx)*tex(i,j,k+1)
c
         tvt2mx=-vdtsmx*areat(i,j,k+1)*2.d0
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        *(tb(i,j,k)-tb(i,j,k+1))
         svt2mx=-vdtsmx*areat(i,j,k+1)*2.d0
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        *(sb(i,j,k)-sb(i,j,k+1))
c
         tvt2=-areat(i,j,k+1)*(
     &        vddt(i,j,k+1)*(tb(i,j,k)-tb(i,j,k+1))*2.d0
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        +( sxw*( aip*(1.-ar)/(1.+sxw*sxw+sys*sys)*0.5
     &        +aip*(1.-ar)/(1.+sxw*sxw+syn*syn)*0.5+agm )
     &        *( (dzt(i,j,k+1)*tb(i,j,k)+dzt(i,j,k)*tb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        -(dzt(i-1,j,k+1)*tb(i-1,j,k)
     &        +dzt(i-1,j,k)*tb(i-1,j,k+1))
     &        /(dzt(i-1,j,k)+dzt(i-1,j,k+1)+1.-tex(i-1,j,k)) )
     &        *(dzumin(i,j-1,k)+dzumin(i,j,k)+dzumin(i,j-1,k+1)
     &        +dzumin(i,j,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        +sxe*( aip*(1.-ar)/(1.+sxe*sxe+sys*sys)*0.5
     &        +aip*(1.-ar)/(1.+sxe*sxe+syn*syn)*0.5+agm )
     &        *( (dzt(i+1,j,k+1)*tb(i+1,j,k)
     &        +dzt(i+1,j,k)*tb(i+1,j,k+1))
     &        /(dzt(i+1,j,k)+dzt(i+1,j,k+1)+1.-tex(i+1,j,k))
     &        -(dzt(i,j,k+1)*tb(i,j,k)+dzt(i,j,k)*tb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k)) )
     &        *(dzumin(i+1,j-1,k)+dzumin(i+1,j,k)+dzumin(i+1,j-1,k+1)
     &        +dzumin(i+1,j,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        )*dxr*cstr(j)
     &        +( sys*( aip*(1.-ar)/(1.+sxw*sxw+sys*sys)*0.5
     &        +aip*(1.-ar)/(1.+sxe*sxe+sys*sys)*0.5+agm )
     &        *( (dzt(i,j,k+1)*tb(i,j,k)+dzt(i,j,k)*tb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        -(dzt(i,j-1,k+1)*tb(i,j-1,k)
     &        +dzt(i,j-1,k)*tb(i,j-1,k+1))
     &        /(dzt(i,j-1,k)+dzt(i,j-1,k+1)+1.-tex(i,j-1,k)) )
     &        *(dzvmin(i-1,j,k)+dzvmin(i,j,k)+dzvmin(i-1,j,k+1)
     &        +dzvmin(i,j,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        +syn*( aip*(1.-ar)/(1.+sxw*sxw+syn*syn)*0.5
     &        +aip*(1.-ar)/(1.+sxe*sxe+syn*syn)*0.5+agm )
     &        *( (dzt(i,j+1,k+1)*tb(i,j+1,k)
     &        +dzt(i,j+1,k)*tb(i,j+1,k+1))
     &        /(dzt(i,j+1,k)+dzt(i,j+1,k+1)+1.-tex(i,j+1,k))
     &        -(dzt(i,j,k+1)*tb(i,j,k)+dzt(i,j,k)*tb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k)) )
     &        *(dzvmin(i-1,j+1,k)+dzvmin(i,j+1,k)+dzvmin(i-1,j+1,k+1)
     &        +dzvmin(i,j+1,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        )*dyr )
c
         svt2=-areat(i,j,k+1)*(
     &        vdds(i,j,k+1)*(sb(i,j,k)-sb(i,j,k+1))*2.d0
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        +( sxw*( aip*(1.-ar)/(1.+sxw*sxw+sys*sys)*0.5
     &        +aip*(1.-ar)/(1.+sxw*sxw+syn*syn)*0.5+agm )
     &        *( (dzt(i,j,k+1)*sb(i,j,k)+dzt(i,j,k)*sb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        -(dzt(i-1,j,k+1)*sb(i-1,j,k)
     &        +dzt(i-1,j,k)*sb(i-1,j,k+1))
     &        /(dzt(i-1,j,k)+dzt(i-1,j,k+1)+1.-tex(i-1,j,k)) )
     &        *(dzumin(i,j-1,k)+dzumin(i,j,k)+dzumin(i,j-1,k+1)
     &        +dzumin(i,j,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        +sxe*( aip*(1.-ar)/(1.+sxe*sxe+sys*sys)*0.5
     &        +aip*(1.-ar)/(1.+sxe*sxe+syn*syn)*0.5+agm )
     &        *( (dzt(i+1,j,k+1)*sb(i+1,j,k)
     &        +dzt(i+1,j,k)*sb(i+1,j,k+1))
     &        /(dzt(i+1,j,k)+dzt(i+1,j,k+1)+1.-tex(i+1,j,k))
     &        -(dzt(i,j,k+1)*sb(i,j,k)+dzt(i,j,k)*sb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k)) )
     &        *(dzumin(i+1,j-1,k)+dzumin(i+1,j,k)+dzumin(i+1,j-1,k+1)
     &        +dzumin(i+1,j,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        )*dxr*cstr(j)
     &        +( sys*( aip*(1.-ar)/(1.+sxw*sxw+sys*sys)*0.5
     &        +aip*(1.-ar)/(1.+sxe*sxe+sys*sys)*0.5+agm )
     &        *( (dzt(i,j,k+1)*sb(i,j,k)+dzt(i,j,k)*sb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k))
     &        -(dzt(i,j-1,k+1)*sb(i,j-1,k)
     &        +dzt(i,j-1,k)*sb(i,j-1,k+1))
     &        /(dzt(i,j-1,k)+dzt(i,j-1,k+1)+1.-tex(i,j-1,k)) )
     &        *(dzvmin(i-1,j,k)+dzvmin(i,j,k)+dzvmin(i-1,j,k+1)
     &        +dzvmin(i,j,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        +syn*( aip*(1.-ar)/(1.+sxw*sxw+syn*syn)*0.5
     &        +aip*(1.-ar)/(1.+sxe*sxe+syn*syn)*0.5+agm )
     &        *((dzt(i,j+1,k+1)*sb(i,j+1,k)
     &        +dzt(i,j+1,k)*sb(i,j+1,k+1))
     &        /(dzt(i,j+1,k)+dzt(i,j+1,k+1)+1.-tex(i,j+1,k))
     &        -(dzt(i,j,k+1)*sb(i,j,k)+dzt(i,j,k)*sb(i,j,k+1))
     &        /(dzt(i,j,k)+dzt(i,j,k+1)+1.-tex(i,j,k)) )
     &        *(dzvmin(i-1,j+1,k)+dzvmin(i,j+1,k)+dzvmin(i-1,j+1,k+1)
     &        +dzvmin(i,j+1,k+1))/(2.*dz(k)+2.*dz(k+1))
     &        )*dyr )
c
         tvt2=dsign(dmin1(dabs(tvt2),dabs(tvt2mx)),tvt2)
         svt2=dsign(dmin1(dabs(svt2),dabs(svt2mx)),svt2)
c
        endif
        endif
        endif
c
        ta(i,j,k)=ta(i,j,k)+tvt2
        ta(i,j,k+1)=ta(i,j,k+1)-tvt2
        sa(i,j,k)=sa(i,j,k)+svt2
        sa(i,j,k+1)=sa(i,j,k+1)-svt2

#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59 .and.
     &   i==33 .and. j==51 ) then
           write(*,*) nkai,'gm-vt',k,ta(i,j,k),sa(i,j,k),tvt2,svt2,
     &     vdts(i,j,k+1)
      endif
#endif

      endif
      enddo
      enddo
      enddo
c
      call wait_t3d_s1p(tvdf,svdf,1,2)
      call wait_t3d_w1p(tudf,sudf,3,4)
c
      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
         ta(i,j,k)=ta(i,j,k)
     $        -tudf(i,j,k)+tudf(i-1,j,k)-tvdf(i,j,k)+tvdf(i,j-1,k)
         sa(i,j,k)=sa(i,j,k)
     $        -sudf(i,j,k)+sudf(i-1,j,k)-svdf(i,j,k)+svdf(i,j-1,k)
      enddo
      enddo
      enddo
      

#ifdef DEBUG1
      if(nkai>482090 .and. nkai < 482120 .and. ip==59) then
          i=33
          j=51
         do k = 1,km
           write(*,*) nkai,'gm-hor',k,ta(i,j,k),
     &    tudf(i,j,k),tudf(i-1,j,k),tvdf(i,j,k),tvdf(i,j-1,k)
           write(*,*) '  svst ',sa(i,j,k),sudf(i,j,k),sudf(i-1,j,k),
     &   svdf(i,j,k),svdf(i,j-1,k)
         enddo
         write(*,*)
      endif
#endif

      return
      end
