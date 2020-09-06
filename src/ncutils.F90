!                     NetCDF Utilities

!  $Id: ncutils.F90 11 2008-12-11 05:08:04Z ishikawa $


      subroutine create_nc_result(ncunit,ncfile,nctitle)
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
         write(*,*)'nf90_create error in create_nc_result : ',nctitle
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

#ifndef CLIMAT
      ncstat=nf90_def_var(ncunit,'nsec',NF90_INT,dimsc,varid)
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : nsec'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'time from begining of the year')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : nsec(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : nsec(unit)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

#endif

      ncstat=nf90_def_var(ncunit,'u',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: u'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','eastward velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : u(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : u(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_var(ncunit,'v',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: v'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','northward velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : v(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : v(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'w',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: w'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'long_name','vertical velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : w(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : w(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

! *****
#ifndef T3OUT

      ncstat=nf90_def_var(ncunit,'t',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: t'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'potential temperature')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : t(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','Celisius')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : t(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_var(ncunit,'s',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: s'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','salinity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : s(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','psu')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : s(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'rhoo',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: rhoo'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'long_name','potential density')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : rhoo(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','g/cm^3')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : rhoo(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'add_offset',addset)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : rhoo(units)'
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



      ncstat=nf90_def_var(ncunit,'tke',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: tke'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'turbulent kinetic enegy')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : tke(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','erg/cm^3')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : tke(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'depml',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: depml'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','mixed layer depth')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : depml(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : depml(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

#if defined(SGOUT) || defined(EXP2006)

      ncstat=nf90_def_var(ncunit,'vdts',NF90_FLOAT,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: vdts'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     &     'vertical mixing coefficient')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : vdts(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','cm^2/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : vdts(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

#endif

#ifdef ICE

      ncstat=nf90_def_var(ncunit,'aice',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: aice'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','ice concentration')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf_put_att_text error : aice(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'hice',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: hice'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','ice height')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : hice(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : hice(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'uice',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: uice'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'eastward ice velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : uice(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : uice(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'vice',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: vice'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'long_name',
     &     'northward ice velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : vice(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : vice(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

#endif


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
      end subroutine create_nc_result


!-------------------- difinition flux data file -------------------------------

      subroutine create_nc_flux(ncunit,ncfile,nctitle)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none
      integer,intent(out) :: ncunit
      character(len=*),intent(in) :: ncfile,nctitle
      integer :: ncstat,dimid,varid,dimsc,dim2d(3)

      if(ip.eq.imaster) then

      ncstat=nf90_create(ncfile,NF90_CLOBBER,ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_create error in create_nc_flux:',nctitle
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,NF90_GLOBAL,'title',nctitle)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_put_att_text error'
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
      dim2d(1)=dimid

      ncstat=nf90_def_dim(ncunit,'lat',jmg,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: lat'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim2d(2)=dimid


      ncstat=nf90_def_dim(ncunit,'time',NF90_UNLIMITED,dimid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_def_dim error: time'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif
      dim2d(3)=dimid
      dimsc=dimid


      ncstat=nf90_def_var(ncunit,'heat',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: heat'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','surface heat flux')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : heat(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','W/cm^2')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : heat(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_def_var(ncunit,'fwater',NF90_FLOAT,dim2d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: heat'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'surface fresh water flux')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : fwater(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif



      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : fwater(units)'
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
      end subroutine create_nc_flux


!----- restoring forward result with NetCDF
#ifdef ASSIM
      subroutine create_nc_forward(ncunit,ncfile,nctitle)
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

      ncstat=nf90_create(ncfile,NF90_CLOBBER,ncunit)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*)'nf90_create error'
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
      ncstat=nf90_def_var(ncunit,'u',NF90_DOUBLE,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: u'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','eastward velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : u(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : u(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_var(ncunit,'v',NF90_DOUBLE,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: v'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','northward velocity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att_text error : v(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','cm/sec')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : v(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_var(ncunit,'t',NF90_DOUBLE,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: t'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'potential temperature')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : t(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','Celisius')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : t(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_var(ncunit,'s',NF90_DOUBLE,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: s'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'long_name','salinity')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : s(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_put_att(ncunit,varid,'units','psu')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : s(units)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif

      ncstat=nf90_def_var(ncunit,'hcl',NF90_DOUBLE,dim2d,varid)
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

      ncstat=nf90_def_var(ncunit,'tke',NF90_DOUBLE,dim3d,varid)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_def_var error: tke'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'long_name',
     $     'turbulent kinetic enegy')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : tke(long_name)'
         write(*,*) nf90_strerror(ncstat)
         stop
      endif
#endif


      ncstat=nf90_put_att(ncunit,varid,'units','erg/cm^3')
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) 'nf90_put_att error : tke(units)'
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
      end subroutine create_nc_forward

#endif
!                    NetCDF I/O interface

!--------------------- 2d data (r8 -> r8) ------------------------

      subroutine write_nc_2d8_r8(arr,ncunit,nt_count,var_name)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      integer,intent(in) ::  ncunit,nt_count
      character(len=*),intent(in) ::  var_name
      real(8),intent(in) :: arr(im,jm)
      integer :: i,j,n,ipxx,ipyy,m,ii,jj,ierr
      integer :: nc_start(3),nc_count(3)
      integer :: ncstat,ncvar
      real*8 x2d8
      common /iocom2d8/ x2d8(img,jmg)
      real(8) :: sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)

!      nc_start(1)=1
!      nc_start(2)=1
!      nc_start(3)=nt_count                       ! time

!      nc_count(1)=img
!      nc_count(2)=jmg
!      nc_count(3)=1                       ! time


      if(ip.eq.imaster) then
        ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) var_name
         write(*,*) nf90_strerror(ncstat)
         return
      endif
#endif

          x2d8(:,:) = 0.
      endif

#ifdef NO_MPI
        do j = 1,jmg
        do i = 1,img
          x2d8(i,j) = arr(i,j)
        enddo
        enddo
        ncstat= nf90_put_var(ncunit,ncvar,x2d4,start=(/ 1,1,nt_count /))
!        ncstat= nf90_put_var(ncunit,ncvar,x2d4,
!     $       start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
         ncstat=nf_close(ncunit)
         stop
      endif
#endif

#else

      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = arr(i,j)
      enddo
      enddo

      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_double_precision,
     &     rbuf,(im-4)*(jm-4),mpi_double_precision,
     &     imaster,comm,ierr)
            if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif

      if(ip.eq.imaster) then
      do m = 1,np
         ipxx = mod(m-1,ipe)
         ipyy = (m-1)/ipe
         do j = 3,jm-2
         do i = 3,im-2
            n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
            ii = i-1 + lw(ipxx)
            jj = j-1 + ls(ipyy)
            if(ii.le.img .and. jj.le.jmg) then
               x2d8(ii,jj) = rbuf(n)
            endif
         enddo
         enddo
      enddo

        ncstat= nf90_put_var(ncunit,ncvar,x2d8,start=(/ 1,1,nt_count /))
#ifdef DEBUG
!      write(*,*) 'var_name : ',var_name
      if(ncstat/=0) then
         write(*,*) var_name,' error :',nf90_strerror(ncstat)
         stop
      endif
#endif
!      write(*,*) 'write hcl exit :', nt_count
      endif
#endif
#endif
      return
      end subroutine write_nc_2d8_r8

!--------------------- 2d data (r8 -> r4) ------------------------

      subroutine write_nc_2d8_r4(arr,ncunit,nt_count,var_name)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      integer,intent(in) ::  ncunit,nt_count
      character(len=*),intent(in) ::  var_name
      real(8),intent(in) :: arr(im,jm)
      integer :: i,j,n,ipxx,ipyy,m,ii,jj,ierr
      integer :: nc_start(3),nc_count(3)
      integer :: ncstat,ncvar
      real(4) :: x2d4
      real(4) :: sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      common /iocom2d4/ x2d4(img,jmg)

!      nc_start(1)=1
!      nc_start(2)=1
!      nc_start(3)=nt_count                       ! time

!      nc_count(1)=img
!      nc_count(2)=jmg
!      nc_count(3)=1                       ! time


      if(ip.eq.imaster) then
        ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) var_name
         write(*,*) nf90_strerror(ncstat)
         return
      endif
#endif

        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = 0.
        enddo
        enddo
      endif

#ifdef NO_MPI
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = arr(i,j)
        enddo
        enddo
        ncstat= nf90_put_var(ncunit,ncvar,x2d4,start=(/ 1,1,nt_count /))
!        ncstat= nf90_put_var(ncunit,ncvar,x2d4,
!     $       start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
         ncstat=nf_close(ncunit)
         stop
      endif
#endif

#else

      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = real(arr(i,j))
      enddo
      enddo

      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_real,
     &     rbuf,(im-4)*(jm-4),mpi_real,imaster,comm,ierr)

      if(ip.eq.imaster) then
      do m = 1,np
         ipxx = mod(m-1,ipe)
         ipyy = (m-1)/ipe
         do j = 3,jm-2
         do i = 3,im-2
            n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
            ii = i-1 + lw(ipxx)
            jj = j-1 + ls(ipyy)
            if(ii.le.img .and. jj.le.jmg) then
               x2d4(ii,jj) = rbuf(n)
            endif
         enddo
         enddo
      enddo

        ncstat= nf90_put_var(ncunit,ncvar,x2d4,start=(/ 1,1,nt_count /))
!      ncstat= nf90_put_var(ncunit,ncvar,x2d4,
!     $     start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
         ncstat=nf90_close(ncunit)
      endif
#endif

      endif
#endif

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
#endif
      return
      end subroutine write_nc_2d8_r4

!--------------------- 2d data (r4) ------------------------

      subroutine write_nc_2d4_r4(arr,ncunit,nt_count,var_name)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      integer,intent(in) ::  ncunit,nt_count
      character(len=*),intent(in) ::  var_name
      integer :: i,j,n,ipxx,ipyy,m,ii,jj,ierr
      integer :: nc_start(3),nc_count(3)
      integer :: ncstat,ncvar
      real(4) :: x2d4
      real(4) :: arr(im,jm)
      real(4) :: sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      common /iocom2d4/ x2d4(img,jmg)

!      nc_start(1)=1
!      nc_start(2)=1
!      nc_start(3)=nt_count                       ! time

!      nc_count(1)=img
!      nc_count(2)=jmg
!      nc_count(3)=1                              ! time

      if(ip.eq.imaster) then
        ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) var_name
         write(*,*) nf90_strerror(ncstat)
         return
      endif
#endif

        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = 0.
        enddo
        enddo
      endif

#ifdef NO_MPI
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = arr(i,j)
        enddo
        enddo

        ncstat= nf90_put_var(ncunit,ncvar,x2d4,start=(/ 1,1,nt_count /))
!        ncstat= nf90_put_var(ncunit,ncvar,x2d4,
!     $       start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
      endif
#endif

#else

      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = real(arr(i,j))
      enddo
      enddo

      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_real,
     &     rbuf,(im-4)*(jm-4),mpi_real,imaster,comm,ierr)

      if(ip.eq.imaster) then
      do m = 1,np
         ipxx = mod(m-1,ipe)
         ipyy = (m-1)/ipe
         do j = 3,jm-2
         do i = 3,im-2
            n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
            ii = i-1 + lw(ipxx)
            jj = j-1 + ls(ipyy)
            if(ii.le.img .and. jj.le.jmg) then
               x2d4(ii,jj) = rbuf(n)
            endif
         enddo
         enddo
      enddo

      ncstat= nf90_put_var(ncunit,ncvar,x2d4,start=(/ 1,1,nt_count /))
!      ncstat= nf90_put_var(ncunit,ncvar,x2d4,
!     $     start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
      endif
#endif

      endif
#endif

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
#endif
      return
      end subroutine write_nc_2d4_r4
!------------------------ 3D data (r8 -> r8) ------------------------
      subroutine write_nc_3d8_r8(arr,ncunit,nt_count,var_name)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      real(8),intent(in) :: arr(im,jm,0:km+1)
      integer,intent(in) :: ncunit,nt_count
      character(len=*),intent(in) ::  var_name
      integer :: i,j,k,n,ipxx,ipyy,m,ii,jj,ierr
      integer :: nc_start(4),nc_count(4)
      integer :: ncstat,ncvar
      real(8) :: sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real*8 x3d8
      common /iocom3d8/ x3d8(img,jmg,0:km+1)

!      nc_start(1)=1
!      nc_start(2)=1
!      nc_start(3)=1
!      nc_start(4)=nt_count

!      nc_count(1)=img
!      nc_count(2)=jmg
!      nc_count(3)=km
!      nc_count(4)=1


!      if(nt_count>28) then
!         write(*,*) ip,nt_count,var_name,'  arr in nc ',maxval(arr)
!      endif
      if(ip==imaster) then
!        write(*,*) 'inq_varid : ',var_name
        ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
!        write(*,*) 'var_id',ncvar
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) var_name,' inq_varid :',nf90_strerror(ncstat)
         return
      endif
#endif
          x3d8(:,:,:) = 0.
      endif

#ifdef NO_MPI
        do j = 1,jmg
        do k = 1,km
        do i = 1,img
          x3d8(i,j,k) = arr(i,j,k)
        enddo
        enddo
        enddo
        ncstat= nf90_put_var(ncunit,ncvar,x3d8(:,:,1:km),
     $       start=(/1,1,1,nt_count /))
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
         ncstat=nf90_close(ncunit)
         stop
      endif
#endif

#else

      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
      enddo

      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_double_precision,
     &     rbuf,(im-4)*(jm-4)*km,mpi_double_precision,
     &     imaster,comm,ierr)
      call mpi_barrier(comm,ierr)
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif

!      write(*,*) 'mpi is OK',ip,nt_count,var_name
      if(ip==imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(k-1)+(im-4)*km*(j-3)
     &        +(im-4)*(jm-4)*km*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
          if(ii.le.img .and. jj.le.jmg) then
            x3d8(ii,jj,k) = rbuf(n)
          endif
        enddo
        enddo
        enddo
      enddo

!      write(*,*) 'var_name : ',var_name,ncvar,nt_count,maxval(x3d8)
      ncstat= nf90_put_var(ncunit,ncvar,x3d8(:,:,1:km),
     &     start=(/1,1,1,nt_count/))
#ifdef DEBUG
!      write(*,*) var_name,' is written:',ncstat
      if(ncstat/=0) then
         write(*,*) var_name,'  putvar: ',nf90_strerror(ncstat)
      endif
#endif

      endif
#endif

#endif
!      write(*,*) 'write nc exit',ip
      return
      end subroutine write_nc_3d8_r8

!------------------------ 3D data (r8 -> r4) ------------------------
      subroutine write_nc_3d8_r4(arr,ncunit,nt_count,var_name)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      real(8),intent(in) :: arr(im,jm,0:km+1)
      integer,intent(in) :: ncunit,nt_count
      character(len=*),intent(in) ::  var_name
      integer :: i,j,k,n,ipxx,ipyy,m,ii,jj,ierr
      integer :: nc_start(4),nc_count(4)
      integer :: ncstat,ncvar
      real(4) :: sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real(4) :: x3d4
      common /iocom3d4/ x3d4(img,jmg,0:km+1)

!      nc_start(1)=1
!      nc_start(2)=1
!      nc_start(3)=1
!      nc_start(4)=nt_count

!      nc_count(1)=img
!      nc_count(2)=jmg
!      nc_count(3)=km
!      nc_count(4)=1


      if(ip.eq.imaster) then
        ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) var_name
         write(*,*) nf90_strerror(ncstat)
         return
      endif
#endif

        do k = 0,km+1
        do j = 1,jmg
        do i = 1,img
          x3d4(i,j,k) = 0.
        enddo
        enddo
        enddo
      endif

#ifdef NO_MPI
        do j = 1,jmg
        do k = 1,km
        do i = 1,img
          x3d4(i,j,k) = arr(i,j,k)
        enddo
        enddo
        enddo
        ncstat= nf90_put_var(ncunit,ncvar,x3d4(:,:,1:km),
     $       start=(/1,1,1,nt_count /))
!        ncstat= nf90_put_var(ncunit,ncvar,x3d4(:,:,1:km),
!     $       start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
         ncstat=nf90_close(ncunit)
         stop
      endif
#endif

#else

      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = real(arr(i,j,k))
      enddo
      enddo
      enddo

      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_real,
     &     rbuf,(im-4)*(jm-4)*km,mpi_real,imaster,comm,ierr)

      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(k-1)+(im-4)*km*(j-3)
     &        +(im-4)*(jm-4)*km*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
          if(ii.le.img .and. jj.le.jmg) then
            x3d4(ii,jj,k) = rbuf(n)
          endif
        enddo
        enddo
        enddo
      enddo

      ncstat= nf90_put_var(ncunit,ncvar,x3d4(:,:,1:km),
     $     start=(/1,1,1,nt_count /))
!      ncstat= nf90_put_var(ncunit,ncvar,x3d4(:,:,1:km),
!     $     start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
      endif
#endif

      endif
#endif

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
#endif
      return
      end subroutine write_nc_3d8_r4

!------------------------ 3D data (r4) --------------------------
      subroutine write_nc_3d4_r4(arr,ncunit,nt_count,var_name)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      integer,intent(in) ::  ncunit,nt_count
      character(len=*),intent(in) ::  var_name
      integer :: nc_start(4),nc_count(4)
      integer :: i,j,k,n,ipxx,ipyy,m,ii,jj,ierr
      integer :: ncstat,ncvar
      real(4) ::  arr(im,jm,0:km+1)
      real(4) ::  sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real(4) ::  x3d4
      common /iocom3d4/ x3d4(img,jmg,0:km+1)

 !     nc_start(1)=1
 !     nc_start(2)=1
 !     nc_start(3)=1
 !     nc_start(4)=nt_count

 !     nc_count(1)=img
 !     nc_count(2)=jmg
 !     nc_count(3)=km
 !     nc_count(4)=1


      if(ip.eq.imaster) then
        ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) var_name
         write(*,*) nf90_strerror(ncstat)
         return
      endif
#endif

        do k = 0,km+1
        do j = 1,jmg
        do i = 1,img
          x3d4(i,j,k) = 0.
        enddo
        enddo
        enddo
      endif
c
#ifdef NO_MPI
        do j = 1,jmg
        do k = 1,km
        do i = 1,img
          x3d4(i,j,k) = arr(i,j,k)
        enddo
        enddo
        enddo
        ncstat= nf90_put_var(ncunit,ncvar,x3d4(:,:,1:km),
     $       start=(/1,1,1,nt_count /))
!       ncstat= nf90_put_var(ncunit,ncvar,x3d4,
!     $       start=nc_start,count=nc_count)
#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
      endif
#endif

#else

      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = real(arr(i,j,k))
      enddo
      enddo
      enddo

      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_real,
     &     rbuf,(im-4)*(jm-4)*km,mpi_real,imaster,comm,ierr)

      if(ip.eq.imaster) then
      do m = 1,np
         ipxx = mod(m-1,ipe)
         ipyy = (m-1)/ipe
         do k = 1,km
         do j = 3,jm-2
         do i = 3,im-2
            n = i-2+(im-4)*(k-1)+(im-4)*km*(j-3)
     &           +(im-4)*(jm-4)*km*(m-1)
            ii = i-1 + lw(ipxx)
            jj = j-1 + ls(ipyy)
            if(ii.le.img .and. jj.le.jmg) then
               x3d4(ii,jj,k) = rbuf(n)
            endif
         enddo
         enddo
         enddo
      enddo

      ncstat= nf90_put_var(ncunit,ncvar,x3d4(:,:,1:km),
     $     start=(/1,1,1,nt_count /))
!      ncstat= nf90_put_var(ncunit,ncvar,x3d4(:,:,1:km),
!     $     start=nc_start,count=nc_count)

#ifdef DEBUG
      if(ncstat.ne.0) then
         write(*,*) nf90_strerror(ncstat)
      endif
#endif

      endif

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif

#endif
#endif
      return
      end subroutine write_nc_3d4_r4

! --------------------------------------------------------------
      subroutine read_nc_fd(arr,ncunit,nds,nde,var_name)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      integer,intent(in) :: ncunit,nds,nde
      character(len=*),intent(in) :: var_name
      real(4) :: arr(im,jm,nds:nde)
      integer :: i,j,n,nc_start(3),nc_count(3),ncstat,ncvar
!      real(4) :: x3d4(img,jmg,nds:nde)
      real(4) :: x2d4(img,jmg)

      nc_start(1)=1
      nc_start(2)=1
!      nc_start(3)=nds

      nc_count(1)=img
      nc_count(2)=jmg
!      nc_count(3)=nde-nds+1
      nc_count(3)=1

      if(ip==imaster) then
         ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) var_name
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
      endif

      do n = nds,nde

        nc_start(3)=n

        if(ip== imaster) then
         ncstat=nf90_get_var(ncunit,ncvar,x2d4,start=nc_start,
     $        count=nc_count)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) var_name
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
        endif

!      call bcast_real(x3d4,img*jmg*(nde-nds+1))
      call bcast_real(x2d4,img*jmg)

!      do n=nds,nde
      do j = 1,jml
      do i = 1,iml
!         arr(i,j,n) = x3d4(lw(ip_x)-1+i,ls(ip_y)-1+j,n)
         arr(i,j,n) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo

      enddo

#endif
      return
      end subroutine read_nc_fd

! --------------------------------------------------------------
      subroutine read_nc_3d_r8(arr,ncunit,var_name,nn)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      integer,intent(in) :: ncunit,nn
      character(len=*),intent(in) :: var_name
      real(8) :: arr(im,jm,0:km+1)
      integer :: i,j,k,nc_start(4),nc_count(4),ncstat,ncvar
      real*8 x3d8
      common /iocom3d8/ x3d8(img,jmg,0:km+1)

      nc_start(1)=1
      nc_start(2)=1
      nc_start(3)=1
      nc_start(4)=nn

      nc_count(1)=img
      nc_count(2)=jmg
      nc_count(3)=km
      nc_count(4)=1
!      write(*,*) ip,var_name,' is reading'
      if(ip==imaster) then
         ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) var_name,', inq_varid ',nf90_strerror(ncstat)
            stop
         endif
#endif
      endif

        if(ip== imaster) then
         ncstat=nf90_get_var(ncunit,ncvar,x3d8(:,:,1:km),
     &        start=nc_start,count=nc_count)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) var_name,'get',nf90_strerror(ncstat)
            stop
         endif
#endif
        endif

      call bcast_dble(x3d8,img*jmg*(km+2))

      arr(:,:,:) = 0.
      do k = 1,km
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k) = x3d8(lw(ip_x)-1+i,ls(ip_y)-1+j,k)
      enddo
      enddo
      enddo

!      write(*,*) ip,var_name,' is finish reading'
#endif
      return
      end subroutine read_nc_3d_r8

! --------------------------------------------------------------
      subroutine read_nc_2d_r8(arr,ncunit,var_name,nn)
#ifdef NCIO
      use param
      use mod_mpi
      use netcdf

      implicit none

      integer,intent(in) :: ncunit,nn
      character(len=*),intent(in) :: var_name
      real(8) :: arr(im,jm)
      integer :: i,j,nc_start(3),nc_count(3),ncstat,ncvar
      real*8 x2d8
      common /iocom2d8/ x2d8(img,jmg)

      nc_start(1)=1
      nc_start(2)=1
      nc_start(3)=nn

      nc_count(1)=img
      nc_count(2)=jmg
      nc_count(3)=1

      if(ip==imaster) then
         ncstat=nf90_inq_varid(ncunit,var_name,ncvar)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) var_name
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
      endif

        if(ip== imaster) then
         ncstat=nf90_get_var(ncunit,ncvar,x2d8,start=nc_start,
     &        count=nc_count)
#ifdef DEBUG
         if(ncstat.ne.0) then
            write(*,*) var_name
            write(*,*) nf90_strerror(ncstat)
            stop
         endif
#endif
        endif

      call bcast_dble(x2d8,img*jmg)

      do j = 1,jml
      do i = 1,iml
         arr(i,j) = x2d8(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo

#endif
      return
      end subroutine read_nc_2d_r8
