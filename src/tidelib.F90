ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE write2d                                            c
c
c  $Id: tidelib.F90 2009-06-15 luu $
c                                                                   c
c     write real*4 tidal data                                       c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE write2d(ncfileid)


      use param
      use mod_mpi

#ifdef NCIO
      use netcdf
#endif

      implicit none

#include "common.h"

      real*8 aday,ayear
      integer iyear,imin,ihour,iday,month
      integer ncfileid
#ifdef NCIO
      integer :: ncvar,ncstat
#endif

#ifdef NCIO
      ncnt_2d = ncnt_2d + 1

      if(ip==imaster) then

         write(*,*) 'write2d: write tide data at step ',nnmats

         ncstat=nf90_inq_varid(ncfileid,'ahour',ncvar)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'ahour'
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

         ncstat=nf90_put_var(ncfileid,ncvar,ahour,ncnt_2d)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

      endif !(ip==imaster)

      call write_nc_2d8_r4(hcl,ncfileid,ncnt_2d,'hcl')
      call write_nc_2d8_r4(sfu,ncfileid,ncnt_2d,'sfu')
      call write_nc_2d8_r4(sfv,ncfileid,ncnt_2d,'sfv')


      if(ip.eq.imaster) then
         ncstat=nf90_sync(ncfileid)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif
      endif
#else !NCIO
      if(ip.eq.imaster) then
        write(ncfileid) ahour
      endif
      call write_2d8_r4(hcl,ncfileid)
      call write_2d8_r4(sfu,ncfileid)
      call write_2d8_r4(sfv,ncfileid)
#endif !NCIO
c
c
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
        write(6,*)nkai,'step : written on tide ncfileid file'
        write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &       ' hour=',ihour,' min=',imin
        write(6,*)' '
      endif

      return
      END SUBROUTINE write2d




#ifdef EQ2OUT        

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE writeq                                            c
c
c  $Id: tidelib.F90 2009-06-15 luu $
c                                                                   c
c     write real*4 tidal data                                       c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE writeq(ncfileid)


      use param
      use mod_mpi

#ifdef NCIO
      use netcdf
#endif

      implicit none

#include "common.h"

      real*8 aday,ayear
      integer iyear,imin,ihour,iday,month
      integer ncfileid
#ifdef NCIO
      integer :: ncvar,ncstat
#endif

#ifdef NCIO
      ncnt_eq = ncnt_eq + 1

      if(ip==imaster) then

         write(*,*) 'writeq: write tide data at step ',nnmats

         ncstat=nf90_inq_varid(ncfileid,'ahour',ncvar)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'ahour'
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

         ncstat=nf90_put_var(ncfileid,ncvar,ahour,ncnt_eq)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

      endif !(ip==imaster)

      call write_nc_2d8_r4(eq_ht,ncfileid,ncnt_eq,'eq_ht')
      call write_nc_2d8_r4(eq_uv,ncfileid,ncnt_eq,'eq_uv')	  
      call write_nc_2d8_r4(eq_id,ncfileid,ncnt_eq,'eq_id')


      if(ip.eq.imaster) then
         ncstat=nf90_sync(ncfileid)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif
      endif
#else !NCIO

#endif !NCIO
c
c
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
        write(6,*)nkai,'step : written on tide ncfileid file'
        write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &       ' hour=',ihour,' min=',imin
        write(6,*)' '
      endif

      return
      END SUBROUTINE writeq

#endif







ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE writeen                                            c
c
c  $Id: tidelib.F90 11 2009-04-15 ishikawa & luu $
c                                                                   c
c     write real*4 energy data                                      c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE writeen


      use param
      use mod_mpi

#ifdef NCIO
      use netcdf
#endif

      implicit none

#include "common.h"

      real*8 aday,ayear
      integer iyear,imin,ihour,iday,month
#ifdef NCIO
      integer :: ncvar,ncstat
#endif

      double precision ebt(im,jm) !,ebc(im,jm),ea(im,jm),ebcl
      double precision efx(im,jm),efy(im,jm),efxy
      integer i,j,k

      ! kinematic barotropic energy (ebt)
      do i=3,iml-2
      do j=3,jml-2
        ebt(i,j)=.5d0*(sfun(i,j)*sfun(i,j)+sfvn(i,j)*sfvn(i,j))
     $    *areauu(j)*(hrr(i,j)+hclu(i,j))
      ebt(i,j)=ebt(i,j)/ttvol
      enddo
      enddo

      ! kinematic baroclinic energy difference from barotropic (ebc)
      ! do i=3,iml-2
      ! do j=3,jml-2
      ! ebcl=0.d0
      ! do k = 1,km
      !    ebcl=ebcl+.5d0*ex(i,j,k)*((u(i,j,k)-sfun(i,j))**2
      ! $        +(v(i,j,k)-sfvn(i,j))**2)*(areauu(j)*dzu(i,j,k))
      ! enddo
      ! ebc(i,j)=ebcl
      ! enddo
      ! enddo

      ! total depth-averaged energy
      ! do i=3,iml-2
      ! do j=3,jml-2
      ! ea(i,j) = 0.d0
      ! do k = 1,km
      !  ea(i,j)=ea(i,j)+.5d0*ex(i,j,k)
      ! &        *(hrr(i,j)+hu(i,j)+(1.d0-ex(i,j,1)))
      ! &        *(u(i,j,k)**2+v(i,j,k)**2 )*dzu(i,j,k)
      ! enddo
      ! ea(i,j) = ea(i,j)/ttvol
      ! enddo
      ! enddo

      ! energy flux
      do i=3,iml-2
      do j=3,jml-2
        efxy=.5d0*(sfun(i,j)*sfun(i,j)+sfvn(i,j)*sfvn(i,j)+
     $    2.d0*grav*hclu(i,j))*(hrr(i,j)+hclu(i,j))*rho_w
      efx(i,j)=efxy*sfun(i,j)
      efy(i,j)=efxy*sfvn(i,j)
      enddo
      enddo


#ifdef NCIO
      ncnt_en = ncnt_en + 1

      if(ip==imaster) then

         write(*,*) 'writeen: write energy data at step ',nnmats

         ncstat=nf90_inq_varid(recfile_en,'ahour',ncvar)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'ahour'
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

         ncstat=nf90_put_var(recfile_en,ncvar,ahour,ncnt_en)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

      endif !(ip==imaster)

      call write_nc_2d8_r4(ebt,recfile_en,ncnt_en,'ebt')
      ! call write_nc_2d8_r4(efx,recfile_en,ncnt_en,'efx')
      ! call write_nc_2d8_r4(efy,recfile_en,ncnt_en,'efy')

      if(ip.eq.imaster) then
         ncstat=nf90_sync(recfile_en)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif
      endif
#else !NCIO
      if(ip.eq.imaster) then
        write(recfile_en) ahour
      endif
      call write_2d8_r4(ebt,recfile_en)
      ! call write_2d8_r4(efx,recfile_en)
      ! call write_2d8_r4(efy,recfile_en)

#endif !NCIO
c
c
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
        write(6,*)nkai,'step : written energy on recfile_en file'
        write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &       ' hour=',ihour,' min=',imin
        write(6,*)' '
      endif

      return
      END SUBROUTINE writeen




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE writecn                                            c
c
c     save for conitinue file (from open to close)                  c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE writecn (conid,conname)

      use param
      use mod_mpi

      implicit none
#include "common.h"

      real*8 aday,ayear
      integer iyear,month,iday,ihour,imin,isec,conid
      character(len=80) :: conname

      if(ip==imaster) then
         write(*,*) 'continue file : ',conname,conid
         open(conid,trim(conname),form='unformatted')
      endif

      ! nkai=nkai-1
      if(ip.eq.imaster) then
         write(conid) nkai,ahour
         write(conid) pd,pm,ddmna,dmn
      endif
      call write_3d8_r8(u,conid)
      call write_3d8_r8(v,conid)
      call write_3d8_r8(t,conid)
      call write_3d8_r8(s,conid)
      call write_3d8_r8(tke,conid)
      call write_2d8_r8(hcl,conid)
      call write_2d8_r8(um,conid)
      call write_2d8_r8(vm,conid)

      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30
      if(ip==imaster) then
         close(conid)
      endif

      if(ip.eq.imaster) then
         write(*,*)nkai,' step : written on conid file'
         write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &        ' hour=',ihour,' min=',imin
         write(*,*)' '
      endif
      return
      END






ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE create_nc_2d                                       c
c
c     luu
c                                                                   c
c     check energy to validate warm-up                              c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE create_nc_2d(ncunit,ncfile,nctitle)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none
      integer,intent(out) :: ncunit
      character(len=*),intent(in) :: ncfile,nctitle
      integer :: ncstat,dimid,varid,dimsc,dim2d(3),dim3d(4)
      real :: addset=1.


!        NetCDF file setting

      if(ip.eq.imaster) then

!------------------  file create ------------------------------

      ncstat=nf90_create(ncfile,NF90_64BIT_OFFSET,ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_create error in create_nc_2d : ',nctitle
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

!------------------  difinition mode --------------------------


      ncstat=nf90_put_att(ncunit,NF90_GLOBAL,'title',nctitle)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_put_att error title'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_dim(ncunit,'lon',img,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lon'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(1)=dimid
      dim2d(1)=dimid

      ncstat=nf90_def_dim(ncunit,'lat',jmg,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lat'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(2)=dimid
      dim2d(2)=dimid

      ncstat=nf90_def_dim(ncunit,'level',km,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_dim error: level'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(3)=dimid

      ncstat=nf90_def_dim(ncunit,'time',NF90_UNLIMITED,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: time'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(4)=dimid
      dim2d(3)=dimid
      dimsc=dimid


!-------  difinition variable

      ncstat=nf90_def_var(ncunit,'ahour',NF90_DOUBLE,dimsc,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: ahour'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','calculation time')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : ahour(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','hour')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : ahour(unit)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'hcl',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: hcl'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','sea surface height')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf_put_att_text error : hcl(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : hcl(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'sfu',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: sfu'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'eastward barotropic velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : sfu(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : sfu(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_var(ncunit,'sfv',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: sfv'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'northward barotropic velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : sfv(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : sfv(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

!   end definition

      ncstat=nf90_enddef(ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_enddef error'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      endif

#endif
      return
      END SUBROUTINE create_nc_2d




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE create_nc_en                                       c
c
c     luu
c                                                                   c
c     check energy to validate warm-up                              c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!  $Id: ncutils.F90 11 2009-06-15 luu $


      SUBROUTINE create_nc_en(ncunit,ncfile,nctitle)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none
      integer,intent(out) :: ncunit
      character(len=*),intent(in) :: ncfile,nctitle
      integer :: ncstat,dimid,varid,dimsc,dim2d(3),dim3d(4)
      real :: addset=1.


!        NetCDF file setting


      if(ip.eq.imaster) then

!------------------  file create ------------------------------

      ncstat=nf90_create(ncfile,NF90_64BIT_OFFSET,ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_create error in create_nc_en : ',nctitle
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

!------------------  difinition mode --------------------------


      ncstat=nf90_put_att(ncunit,NF90_GLOBAL,'title',nctitle)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_put_att error title'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_dim(ncunit,'lon',img,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lon'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(1)=dimid
      dim2d(1)=dimid

      ncstat=nf90_def_dim(ncunit,'lat',jmg,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lat'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(2)=dimid
      dim2d(2)=dimid

      ncstat=nf90_def_dim(ncunit,'level',km,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_dim error: level'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(3)=dimid

      ncstat=nf90_def_dim(ncunit,'time',NF90_UNLIMITED,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: time'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(4)=dimid
      dim2d(3)=dimid
      dimsc=dimid


!-------  difinition variable

      ncstat=nf90_def_var(ncunit,'ahour',NF90_DOUBLE,dimsc,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: ahour'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','calculation time')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : ahour(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','hour')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : ahour(unit)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'ebt',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: ebt'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','barotropic energy')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf_put_att_text error : ebt(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','g*cm^2/s^2')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : ebt(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


!      ncstat=nf90_def_var(ncunit,'efx',NF90_FLOAT,dim2d,varid)
! #ifdef DEBUG
!       if(ncstat.ne.0) then
!          write(*,*) 'nf90_def_var error: efx'
!          write(*,*) nf90_strerror(ncstat)
!          stop
!      endif
! #endif

!      ncstat=nf90_put_att(ncunit,varid,'long_name',
!     $     'total energy')
!#ifdef DEBUG
!      if(ncstat.ne.0) then
!         write(*,*) 'nf90_put_att error : efx(long_name)'
!         write(*,*) nf90_strerror(ncstat)
!         stop
!      endif
!#endif

!      ncstat=nf90_put_att(ncunit,varid,'units','g/s^3')
!#ifdef DEBUG
!      if(ncstat.ne.0) then
!         write(*,*) 'nf90_put_att error : efx(units)'
!         write(*,*) nf90_strerror(ncstat)
!         stop
!      endif
!#endif


!      ncstat=nf90_def_var(ncunit,'efy',NF90_FLOAT,dim2d,varid)
!#ifdef DEBUG
!      if(ncstat.ne.0) then
!         write(*,*) 'nf90_def_var error: efy'
!         write(*,*) nf90_strerror(ncstat)
!         stop
!      endif
!#endif

!      ncstat=nf90_put_att(ncunit,varid,'long_name',
!     $     'total energy')
!#ifdef DEBUG
!      if(ncstat.ne.0) then
!         write(*,*) 'nf90_put_att error : efy(long_name)'
!         write(*,*) nf90_strerror(ncstat)
!         stop
!      endif
!#endif

!      ncstat=nf90_put_att(ncunit,varid,'units','g/s^3')
!#ifdef DEBUG
!      if(ncstat.ne.0) then
!         write(*,*) 'nf90_put_att error : efy(units)'
!         write(*,*) nf90_strerror(ncstat)
!         stop
!      endif
!#endif


!   end definition

      ncstat=nf90_enddef(ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_enddef error'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      endif

#endif
      return
      END SUBROUTINE create_nc_en




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE create_nc_fast                                     c
c
c     luu
c                                                                   c
c     create fast netcdf file                                       c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!  $Id: ncutils.F90 11 2009-06-15 luu $


      SUBROUTINE create_nc_fast(path,varname,varvalue)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      character(len=80),intent(in) :: path,varname
      character(len=80) :: ncfile
#ifndef BNDTIDE_SAVENAO
      real*8,intent(in) :: varvalue(im,jm)
#else
      real*8,intent(in) :: varvalue(im,jm,32)
#endif
      real :: addset=1.
      integer :: ncstat,dimid,varid,dimsc,dim2d(3),dim3d(4)
      integer :: ncnt_rs_fast,ncunit !must be different with defines in set.F90

!        NetCDF file setting
      if(ip.eq.imaster) then
        ncfile = trim(path)//'rs_bndtide_'//trim(varname)//'.nc'
        write(*,*) 'ncfile',ncfile
        ncunit = 97
!------------------  file create ------------------------------

      ncstat=nf90_create(ncfile,NF90_64BIT_OFFSET,ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_create error in create_nc_fast : ',varname
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

!------------------  difinition mode --------------------------


      ncstat=nf90_def_dim(ncunit,'lon',img,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lon'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(1)=dimid
      dim2d(1)=dimid

      ncstat=nf90_def_dim(ncunit,'lat',jmg,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lat'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(2)=dimid
      dim2d(2)=dimid

      ncstat=nf90_def_dim(ncunit,'level',km,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_dim error: level'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(3)=dimid

      ncstat=nf90_def_dim(ncunit,'time',NF90_UNLIMITED,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: time'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(4)=dimid
      dim2d(3)=dimid
      dimsc=dimid

!-------  difinition variable

#ifndef BNDTIDE_SAVENAO
      ncstat=nf90_def_var(ncunit,trim(varname),NF90_FLOAT,dim2d,varid)
#else
      ncstat=nf90_def_var(ncunit,trim(varname),NF90_FLOAT,dim3d,varid)
#endif

#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: ',varname
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'fast variable')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : (longname) ',varname
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','undefined')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : (units) ',varname
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

!   end definition

      ncstat=nf90_enddef(ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_enddef error'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      endif

#endif


!-------  write variable

      ncnt_rs_fast = 1

#ifdef NCIO

#ifndef BNDTIDE_SAVENAO
      call write_nc_2d8_r4(varvalue,ncunit,ncnt_rs_fast,varname)
#else
      call write_nc_3d8_r4(varvalue,ncunit,ncnt_rs_fast,varname)
#endif

      if(ip.eq.imaster) then
         ncstat=nf90_sync(ncunit)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif
      endif
#else !NCIO

#ifndef BNDTIDE_SAVENAO
      call write_2d8_r4(varvalue,ncunit)
#else
      call write_3d8_r4(varvalue,ncunit)
#endif

#endif !NCIO

      if(ip.eq.imaster) then
         ncstat=nf90_close(ncunit)
      endif


      return
      END SUBROUTINE create_nc_fast




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE savenao                                            c
c
c     luu
c                                                                   c
c     save nao boundary to files                                    c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE savenao(path)

! Produce heights of tidal wave
      use mod_mpi
      use mod_nao
      implicit none

      real*8:: ax(im,jm,tcm),px(im,jm,tcm),nao(im,jm,32),fx(tcm)
      character(len=80) :: path,fname
      integer i,j,k


      call extracttide(ax,px,fx)

      fname = 'nao'

      do i=1,iml
      do j=1,jml
      do k=1,tcm
         nao(i,j,k) = ax(i,j,k)
         nao(i,j,tcm+k) = px(i,j,k)
      enddo
      enddo
      enddo

      call create_nc_fast(path,fname,nao)

      if(ip.eq.imaster)then
         write(*,*)'tidal frequency '
         do k=1,tcm
            write(*,*)k,'f(',k,') = ',fx(k),' cycles/day'
         enddo
      endif

      return
      END SUBROUTINE savenao





#ifdef IDEBUG   

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE create_nc_eqterms                                  c
c
c     luu
c                                                                   c
c     check energy to validate warm-up                              c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE create_nc_eqterms (ncunit,ncfile,nctitle)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none
      integer,intent(out) :: ncunit
      character(len=*),intent(in) :: ncfile,nctitle
      integer :: ncstat,dimid,varid,dimsc,dim2d(3),dim3d(4)
      real :: addset=1.


!        NetCDF file setting

      if(ip.eq.imaster) then

!------------------  file create ------------------------------

      ncstat=nf90_create(ncfile,NF90_64BIT_OFFSET,ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_create error in create_nc_eqterms : ',nctitle
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

!------------------  difinition mode --------------------------


      ncstat=nf90_put_att(ncunit,NF90_GLOBAL,'title',nctitle)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_put_att error title'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_dim(ncunit,'lon',img,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lon'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(1)=dimid
      dim2d(1)=dimid

      ncstat=nf90_def_dim(ncunit,'lat',jmg,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lat'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(2)=dimid
      dim2d(2)=dimid

      ncstat=nf90_def_dim(ncunit,'level',km,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_dim error: level'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(3)=dimid

      ncstat=nf90_def_dim(ncunit,'time',NF90_UNLIMITED,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: time'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim3d(4)=dimid
      dim2d(3)=dimid
      dimsc=dimid


!-------  difinition variable

      ncstat=nf90_def_var(ncunit,'ahour',NF90_DOUBLE,dimsc,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: ahour'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','calculation time')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : ahour(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','hour')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : ahour(unit)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif




      ncstat=nf90_def_var(ncunit,'eq_ht',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: eq_ht'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','eq_ht')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf_put_att_text error : eq_ht(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','unknown')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : eq_ht(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif




      ncstat=nf90_def_var(ncunit,'eq_uv',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: eq_uv'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','eq_uv')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf_put_att_text error : eq_uv(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','unknown')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : eq_uv(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif



      ncstat=nf90_def_var(ncunit,'eq_id',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: eq_id'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','eq_id')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf_put_att_text error : eq_id(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','unknown')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : eq_id(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif



!   end definition

      ncstat=nf90_enddef(ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_enddef error'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      endif

#endif
      return
      END SUBROUTINE create_nc_eqterms





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE writeq                                            c
c
c  $Id: tidelib.F90 2009-06-15 luu $
c                                                                   c
c     write real*4 tidal data                                       c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE writeq(ncfileid)


      use param
      use mod_mpi

#ifdef NCIO
      use netcdf
#endif

      implicit none

#include "common.h"

      real*8 aday,ayear
      integer iyear,imin,ihour,iday,month
      integer ncfileid
#ifdef NCIO
      integer :: ncvar,ncstat
#endif

#ifdef NCIO
      ncnt_eq = ncnt_eq + 1

      if(ip==imaster) then

         write(*,*) 'writeq: write tide data at step ',nnmats

         ncstat=nf90_inq_varid(ncfileid,'ahour',ncvar)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) 'ahour'
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

         ncstat=nf90_put_var(ncfileid,ncvar,ahour,ncnt_eq)

#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif

      endif !(ip==imaster)


      call write_nc_2d8_r4(eq_ht,ncfileid,ncnt_eq,'eq_ht')
      call write_nc_2d8_r4(eq_uv,ncfileid,ncnt_eq,'eq_uv')	  
      call write_nc_2d8_r4(eq_id,ncfileid,ncnt_eq,'eq_id')


      if(ip.eq.imaster) then
         ncstat=nf90_sync(ncfileid)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) nf90_strerror(ncstat)
         endif
#endif
      endif
#endif !NCIO
c
c
      aday=ahour/24.d0
      ayear=aday/dble(nday)
      iyear=ayear+.0000001d0
      imin=(ahour-dble(iyear*nday*24))*60.d0+.01d0
      ihour=imin/60
      imin=imin-ihour*60
      iday=ihour/24
      ihour=ihour-iday*24
      month=iday/30
      iday=iday-month*30

      if(ip.eq.imaster) then
        write(6,*)nkai,'step : written on tide ncfileid file'
        write(*,*)'year=',iyear,' month=',month,' day=',iday,
     &       ' hour=',ihour,' min=',imin
        write(6,*)' '
      endif

      return
      END SUBROUTINE writeq





#endif





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c     SUBROUTINE i2s                                                c
c
c     luu
c                                                                   c
c     convert integer to string                                     c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE i2s (number,maxorder,string)

      implicit none
      integer,intent(in) :: number,maxorder
      character(len=maxorder),intent(out) :: string
      character(len=maxorder) :: temp
      character(1) :: numchar
      integer:: i,n,k

      n = number
      k = 0
      do while (n.gt.0)
         i = mod(n,10)
         n = int(n/10)
         numchar = char(i+48)
         temp = numchar//temp
         k = k + 1
      enddo

      numchar = '0'
      if(k.lt.maxorder)then
         do i = 1,maxorder-k
            temp = numchar//temp
         enddo
      endif
      string = temp

      END