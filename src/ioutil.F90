c
c     i/o utilties for mpi parallelization
c     $Id: ioutil.F90 11 2008-12-11 05:08:04Z ishikawa $
c----------------------------------------------------------
c
      subroutine read_2d_i(iarr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,iarr(im,jm),ibuf_in(img,jmg),i,j
c
      if(ip.eq.imaster) then
        read(nunit) ibuf_in
      endif
c
      call bcast_int(ibuf_in,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        ibuf_in(    1,j) = ibuf_in(img-3,j)
        ibuf_in(    2,j) = ibuf_in(img-2,j)
        ibuf_in(img-1,j) = ibuf_in(    3,j)
        ibuf_in(img  ,j) = ibuf_in(    4,j)
      enddo
#endif
c
      do j = 1,jml
      do i = 1,iml
         iarr(i,j) = ibuf_in(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      return
      end
c
c----------------------------------------------------------
#ifdef GL11M
c
      subroutine read_ref(arr,nunit)

      use param
      use mod_mpi
      implicit none
      integer nunit,i,j,k,n
      real*8 arr(im,jm,0:km+1,0:imn+1)
      real*8 buf_in(img,jmg,0:km+1,0:imn+1)
c
      if(ip.eq.imaster) then
        read(nunit) buf_in
      endif
c
      call bcast_dble(buf_in,img*jmg*(km+2)*(imn+2))
c
#ifdef CYCLIC
      do n = 0,imn+1
      do k = 0,km+1
      do j = 1,jmg
        buf_in(    1,j,k,n) = buf_in(img-3,j,k,n)
        buf_in(    2,j,k,n) = buf_in(img-2,j,k,n)
        buf_in(img-1,j,k,n) = buf_in(    3,j,k,n)
        buf_in(img  ,j,k,n) = buf_in(    4,j,k,n)
      enddo
      enddo
      enddo
#endif
c
      do m = 0,imn+1
      do k = 0,km+1
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k,m) = buf_in(lw(ip_x)-1+i,ls(ip_y)-1+j,k,m)
      enddo
      enddo
      enddo
      enddo
c
      return
      end
c
#endif
c----------------------------------------------------------
c
      subroutine read_ref1(arr,nunit)

      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,k
      real*4 arr(im,jm,0:km+1)
      real*4 x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        read(nunit) x2d4
      endif
c
      call bcast_real(x2d4,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        x2d4(    1,j) = x2d4(img-3,j)
        x2d4(    2,j) = x2d4(img-2,j)
        x2d4(img-1,j) = x2d4(    3,j)
        x2d4(img  ,j) = x2d4(    4,j)
      enddo
#endif
c
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      enddo
c
      do j=1,jml
      do i=1,iml
         arr(i,j,0)=0.
         arr(i,j,km+1)=0.
      enddo
      enddo
c
c
      return
      end
c
c----------------------------------------------------------
#ifdef GL11M
c
      subroutine read_fd(arr,nunit)

      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,n
      real arr(im,jm,nsf)
      real buf_in(img,jmg,nsf)
c
      if(ip.eq.imaster) then
        read(nunit) buf_in
      endif
c
      call bcast_real(buf_in,img*jmg*nsf)
c
#ifdef CYCLIC
      do k = 1,nsf
      do j = 1,jmg
        buf_in(    1,j,k) = buf_in(img-3,j,k)
        buf_in(    2,j,k) = buf_in(img-2,j,k)
        buf_in(img-1,j,k) = buf_in(    3,j,k)
        buf_in(img  ,j,k) = buf_in(    4,j,k)
      enddo
      enddo
#endif
c
      do n = 1,nsf
      do j = 1,jml
      do i = 1,iml
         arr(i,j,n) = buf_in(lw(ip_x)-1+i,ls(ip_y)-1+j,n)
      enddo
      enddo
      enddo
c
c
      return
      end
c
#endif
c----------------------------------------------------------
#ifdef PC68M
c
      subroutine read_fd(arr,nunit)

      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,n
      real arr(im,jm,nsf)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do n=1,nsf
c
      if(ip.eq.imaster) then
        read(nunit) x2d4
      endif
c
      call bcast_real(x2d4,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        x2d4(    1,j) = x2d4(img-3,j)
        x2d4(    2,j) = x2d4(img-2,j)
        x2d4(img-1,j) = x2d4(    3,j)
        x2d4(img  ,j) = x2d4(    4,j)
      enddo
#endif
      do j = 1,jml
      do i = 1,iml
         arr(i,j,n) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      enddo
c
c#ifdef DEBUG
c      jbuf = ls(ip)-1
c      write(*,*) 'forc :',ip,jbuf,arr(151,11,1),arr(153,13,1),
c     &    buf_in(151,jbuf+11,1),buf_in(153,jbuf+13,1)
c#endif
c
      return
      end
c
#endif
c-------------------------------------------------------------
#ifdef NWNPAC
c
      subroutine read_fd(arr,nunit)

      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,n
#ifdef CLIMAT
      real arr(im,jm,nday)
#else
      real arr(im,jm,0:nday+1)
#endif
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do n=1,nday
c
      if(ip.eq.imaster) then
        read(nunit) x2d4
      endif
c
      call bcast_real(x2d4,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        x2d4(    1,j) = x2d4(img-3,j)
        x2d4(    2,j) = x2d4(img-2,j)
        x2d4(img-1,j) = x2d4(    3,j)
        x2d4(img  ,j) = x2d4(    4,j)
      enddo
#endif
      do j = 1,jml
      do i = 1,iml
         arr(i,j,n) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      enddo

      return
      end
c
#endif
c----------------------------------------------------------
#ifdef JP68M
c
      subroutine read_fd(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,n
#ifdef CLIMAT
      real arr(im,jm,nday)
#else
      real arr(im,jm,0:nday+1)
#endif
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do n=1,nday
c
      if(ip.eq.imaster) then
        read(nunit) x2d4
      endif
c
      call bcast_real(x2d4,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        x2d4(    1,j) = x2d4(img-3,j)
        x2d4(    2,j) = x2d4(img-2,j)
        x2d4(img-1,j) = x2d4(    3,j)
        x2d4(img  ,j) = x2d4(    4,j)
      enddo
#endif
      do j = 1,jml
      do i = 1,iml
         arr(i,j,n) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      enddo

      return
      end
c
#endif
c-----------------------------------------------------------------
      subroutine read_fd_1d(arr,nunit)
      use param
      use mod_mpi
      implicit none
c
      integer nunit,i,j
      real arr(im,jm)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      if(ip.eq.imaster) then
        read(nunit) x2d4
      endif
c
      call bcast_real(x2d4,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        x2d4(    1,j) = x2d4(img-3,j)
        x2d4(    2,j) = x2d4(img-2,j)
        x2d4(img-1,j) = x2d4(    3,j)
        x2d4(img  ,j) = x2d4(    4,j)
      enddo
#endif
      do j = 1,jml
      do i = 1,iml
         arr(i,j) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------
#ifdef PC68M
c
      subroutine read_3d_r4(arr,nunit)
	  
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,k
      real arr(im,jm,0:km+1)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        read(nunit) x2d4
      endif
c
      call bcast_real(x2d4,img*jmg)
c
#ifdef CYCLIC
c      do k = 0,km+1
      do j = 1,jmg
        x2d4(    1,j) = x2d4(img-3,j)
        x2d4(    2,j) = x2d4(img-2,j)
        x2d4(img-1,j) = x2d4(    3,j)
        x2d4(img  ,j) = x2d4(    4,j)
      enddo
c      enddo
#endif
c
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      enddo
c
      do j=1,jm
      do i=1,im
         arr(i,j,0)=0.
         arr(i,j,km+1)=0.
      enddo
      enddo
c
      return
      end
c
#endif
c---------------------------------------------------------------------------
#ifdef JP68M
c
      subroutine read_3d_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,k
      real arr(im,jm,0:km+1)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        read(nunit) x2d4
      endif
c
      call bcast_real(x2d4,img*jmg)
c
#ifdef CYCLIC
c      do k = 0,km+1
      do j = 1,jmg
        x2d4(    1,j) = x2d4(img-3,j)
        x2d4(    2,j) = x2d4(img-2,j)
        x2d4(img-1,j) = x2d4(    3,j)
        x2d4(img  ,j) = x2d4(    4,j)
      enddo
c      enddo
#endif
c
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k) = x2d4(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      enddo
c
      do j=1,jml
      do i=1,iml
         arr(i,j,0)=0.
         arr(i,j,km+1)=0.
      enddo
      enddo
c
      return
      end
c
#endif
c----------------------------------------------------------
#if defined(GL11M) || defined(NWNPAC)
c
      subroutine read_3d_r48(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,k
      real*8 arr(im,jm,0:km+1)
      real buf_in(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
        read(nunit) buf_in
      endif
c
      call bcast_real(buf_in,img*jmg*(km+2))
c
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jmg
        buf_in(    1,j,k) = buf_in(img-3,j,k)
        buf_in(    2,j,k) = buf_in(img-2,j,k)
        buf_in(img-1,j,k) = buf_in(    3,j,k)
        buf_in(img  ,j,k) = buf_in(    4,j,k)
      enddo
      enddo
#endif
c
      do k = 0,km+1
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k) = dble( buf_in(lw(ip_x)-1+i,ls(ip_y)-1+j,k) )
      enddo
      enddo
      enddo
c
      return
      end
c
#endif
c----------------------------------------------------------
c
      subroutine read_3d_r8(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,k
      real*8 arr(im,jm,0:km+1)
      real*8 x3d8
      common /iocom3d8/ x3d8(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
        read(nunit) x3d8
      endif
c
      call bcast_dble(x3d8,img*jmg*(km+2))
c
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jmg
        x3d8(    1,j,k) = x3d8(img-3,j,k)
        x3d8(    2,j,k) = x3d8(img-2,j,k)
        x3d8(img-1,j,k) = x3d8(    3,j,k)
        x3d8(img  ,j,k) = x3d8(    4,j,k)
      enddo
      enddo
#endif
c
      do k = 0,km+1
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k) = x3d8(lw(ip_x)-1+i,ls(ip_y)-1+j,k)
      enddo
      enddo
      enddo
c
      return
      end
c
c----------------------------------------------------------
#ifdef PC68M
c
      subroutine read_3d_2d_r8(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j,k
      real*8 arr(im,jm,0:km+1)
      real*8 x2d8
      common /iocom2d8/ x2d8(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        read(nunit) x2d8
      endif
c
      call bcast_dble(x2d8,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        x2d8(    1,j) = x2d8(img-3,j)
        x2d8(    2,j) = x2d8(img-2,j)
        x2d8(img-1,j) = x2d8(    3,j)
        x2d8(img  ,j) = x2d8(    4,j)
      enddo
#endif
c
      do j = 1,jml
      do i = 1,iml
         arr(i,j,k) = x2d8(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      enddo
c
      do j=1,jm
      do i=1,im
         arr(i,j,0)=0.
         arr(i,j,km+1)=0.
      enddo
      enddo
c
      return
      end
c
#endif
c----------------------------------------------------------
c
      subroutine read_2d_r8(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,i,j
      real*8 arr(im,jm)
      real*8 x2d8
      common /iocom2d8/ x2d8(img,jmg)
c
      if(ip.eq.imaster) then
        read(nunit) x2d8
      endif
c
      call bcast_dble(x2d8,img*jmg)
c
#ifdef CYCLIC
      do j = 1,jmg
        x2d8(    1,j) = x2d8(img-3,j)
        x2d8(    2,j) = x2d8(img-2,j)
        x2d8(img-1,j) = x2d8(    3,j)
        x2d8(img  ,j) = x2d8(    4,j)
      enddo
#endif
c
      do j = 1,jml
      do i = 1,iml
         arr(i,j) = x2d8(lw(ip_x)-1+i,ls(ip_y)-1+j)
      enddo
      enddo
c
      return
      end
c
c------------------------------------------------------------------------
c
      subroutine write_3d8_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,k,m,ii,jj,ipxx,ipyy
      real*8 arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real x3d4
      common /iocom3d4/ x3d4(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
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
      write(nunit) x3d4
#else
c
      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = real(arr(i,j,k))
      enddo
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_real,
     &     rbuf,(im-4)*(jm-4)*km,mpi_real,imaster,comm,ierr)
c
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
c
        write(nunit) x3d4
      endif
#endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      return
      end
c
c------------------------------------------------------------------------
#ifdef PC68M
c
      subroutine write_3d4_2d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,k,m,ii,jj,ipxx,ipyy
      real arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = 0.
        enddo
        enddo
      endif
c
#ifdef NO_MPI
        do j = 3,jmg-2
        do i = 3,img-2
          x2d4(i,j) = arr(i,j,k)
        enddo
        enddo
c        enddo
      write(nunit) x2d4
#else
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_real,
     &     rbuf,(im-4)*(jm-4),mpi_real,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
c        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
          if(ii.le.img .and. jj.le.jmg) then
c            x2d4(ii,jj,k) = rbuf(n)
            x2d4(ii,jj) = rbuf(n)
          endif
        enddo
        enddo
c        enddo
      enddo
        write(nunit) x2d4
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      call mpi_barrier(comm,ierr)
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
c
      enddo
c
      return
      end
c
#endif
c------------------------------------------------------------------------
c
      subroutine write_2d8_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,ii,jj,m,k,ipxx,ipyy
      real*8 arr(im,jm)
      real sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = 0.
        enddo
        enddo
      endif
c
#ifdef NO_MPI
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = arr(i,j)
        enddo
        enddo
      write(nunit) x2d4
#else
c
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = arr(i,j)
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_real,
     &     rbuf,(im-4)*(jm-4),mpi_real,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
              ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
!        do k = 1,km
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
!        enddo
      enddo
        write(nunit) x2d4
      endif
#endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      return
      end
c
c
c------------------------------------------------------------------------
c
      subroutine write_3d8_r8(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,ii,jj,k,m,ipxx,ipyy
      real*8 arr(im,jm,0:km+1)
      real*8 sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real*8 x3d8
      common /iocom3d8/ x3d8(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
        do k = 0,km+1
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
      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_double_precision,
     &     rbuf,(im-4)*(jm-4)*km,mpi_double_precision,
     &     imaster,comm,ierr)
c
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
            x3d8(ii,jj,k) = rbuf(n)
          endif
        enddo
        enddo
        enddo
      enddo
        write(nunit) x3d8
      endif
#endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      return
      end
c------------------------------------------------------------------------
c
      subroutine write_2d8_r8(arr,nunit)
      use param
      use mod_mpi
      implicit none
	  
      integer ierr,nunit,n,i,j,n,m,ii,jj,ipxx,ipyy
      real*8 arr(im,jm)
      real*8 sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real*8 x2d8
      common /iocom2d8/ x2d8(img,jmg)
c
      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d8(i,j) = 0.
        enddo
        enddo
      endif
c
#ifdef NO_MPI
      write(nunit) arr
#else
c
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = arr(i,j)
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_double_precision,
     &     rbuf,(im-4)*(jm-4),mpi_double_precision,
     &     imaster,comm,ierr)
c
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
        write(nunit) x2d8
      endif
#endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      return
      end
c
c------------------------------------------------------------------------
#ifdef GL11M
c
      subroutine write_3d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,k,m,ii,jj,ipxx,ipyy
      real arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real x3d4
      common /iocom3d4/ x3d4(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
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
      write(nunit) arr
#else
c
      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_real,
     &     rbuf,(im-4)*(jm-4)*km,mpi_real,imaster,comm,ierr)
c
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
        write(nunit) x3d4
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
#endif
c---------------------------------------------------------------------------
#ifdef JP44
c
      subroutine write_3d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none
	  
      integer ierr,nunit,n,i,j,k,ipxx,ipyy,m,ii,jj
      real arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real x3d4
      common /iocom3d4/ x3d4(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
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
      write(nunit) arr
#else
c
      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_real,
     &     rbuf,(im-4)*(jm-4)*km,mpi_real,imaster,comm,ierr)
c
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
        write(nunit) x3d4
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
#endif
c-----------------------------------------------------------------------------
#ifdef NWNPAC
c
      subroutine write_3d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none
	  
      integer ierr,nunit,n,i,j,k,ipxx,ipyy,m,ii,jj
      real arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real x3d4
      common /iocom3d4/ x3d4(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
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
      write(nunit) arr
#else
c
      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_real,
     &     rbuf,(im-4)*(jm-4)*km,mpi_real,imaster,comm,ierr)
c
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
        write(nunit) x3d4
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
#endif
c-----------------------------------------------------------------------------
#ifdef JP68M
c
      subroutine write_3d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,k,ipxx,ipyy,m,ii,jj
      real arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)*km),rbuf((im-4)*(jm-4)*km*np)
      real x3d4
      common /iocom3d4/ x3d4(img,jmg,0:km+1)
c
      if(ip.eq.imaster) then
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
      write(nunit) arr
#else
c
      do k = 1,km
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(k-1) + (im-4)*km*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4)*km,mpi_real,
     &     rbuf,(im-4)*(jm-4)*km,mpi_real,imaster,comm,ierr)
c
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
        write(nunit) x3d4
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
#endif
c------------------------------------------------------------------------
c
      subroutine write_2d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,m,ii,jj,ipxx,ipyy
      real arr(im,jm)
      real sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = 0.
        enddo
        enddo
      endif
c
#ifdef NO_MPI
      write(nunit) arr
#else
c
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = arr(i,j)
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_real,
     &     rbuf,(im-4)*(jm-4),mpi_real,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
c        do k = 1,km
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
c        enddo
      enddo
        write(nunit) x2d4
      endif
c
#endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      return
      end
c
c------------------------------------------------------------------------
#ifdef PC68M
c
      subroutine write_3d8_2d8_r8(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,m,k,ii,jj,ipxx,ipyy
      real*8 arr(im,jm,0:km+1)
      real*8 sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real*8 x2d8
      common /iocom2d8/ x2d8(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d8(i,j) = 0.
        enddo
        enddo
      endif
c
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_double_precision,
     &     rbuf,(im-4)*(jm-4),mpi_double_precision,
     &     imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
c        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
          if(ii.le.img .and. jj.le.jmg) then
c            x2d8(ii,jj,k) = rbuf(n)
            x2d8(ii,jj) = rbuf(n)
          endif
        enddo
        enddo
c        enddo
      enddo
        write(nunit)x2d8
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      call mpi_barrier(comm,ierr)
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      enddo
c
      return
      end
c
#endif
c--------------------------------------------------------------------
#ifdef NWNPAC
c
      subroutine write_3d8_2d8_r8(arr,nunit)
      use param
      use mod_mpi
      implicit none
c
      integer ierr,nunit,n,i,j,ipxx,ipyy,m,ii,jj,k
      real*8 arr(im,jm,0:km+1)
      real*8 sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real*8 x2d8
      common /iocom2d8/ x2d8(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d8(i,j) = 0.
        enddo
        enddo
      endif
c
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = arr(i,j,k)
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_double_precision,
     &     rbuf,(im-4)*(jm-4),mpi_double_precision,
     &     imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
!        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
          if(ii.le.img .and. jj.le.jmg) then
            x2d8(ii,jj) = rbuf(n)
!            x2d8(ii,jj,k) = rbuf(n)
          endif
        enddo
        enddo
!        enddo
      enddo
        write(nunit)x2d8
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      call mpi_barrier(comm,ierr)
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      enddo
c
      return
      end
c
#endif
c------------------------------------------------------------------------
#ifdef PC68M
c
      subroutine write_3d8_2d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none

      integer ierr,nunit,n,i,j,k,m,ii,jj,ipxx,ipyy
      real*8 arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)

      do k=1,km

      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = 0.
        enddo
        enddo
      endif

      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = real(arr(i,j,k))
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_real,
     &     rbuf,(im-4)*(jm-4),mpi_real,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
c        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
          if(ii.le.img .and. jj.le.jmg) then
c            x2d4(ii,jj,k) = rbuf(n)
            x2d4(ii,jj) = rbuf(n)
          endif
        enddo
        enddo
c        enddo
      enddo
        write(nunit) x2d4
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      ierr=0
      call mpi_barrier(comm,ierr)
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      enddo
c
c
      return
      end
c
#endif
c------------------------------------------------------------------------
#ifdef NWNPAC
c
      subroutine write_3d8_2d4_r4(arr,nunit)
      use param
      use mod_mpi
      implicit none
c
      integer ierr,nunit,n,i,j,ipxx,ipyy,m,ii,jj,k
      real*8 arr(im,jm,0:km+1)
      real sbuf((im-4)*(jm-4)),rbuf((im-4)*(jm-4)*np)
      real x2d4
      common /iocom2d4/ x2d4(img,jmg)
c
      do k=1,km
c
      if(ip.eq.imaster) then
        do j = 1,jmg
        do i = 1,img
          x2d4(i,j) = 0.
        enddo
        enddo
      endif
c
      do j = 3,jm-2
      do i = 3,im-2
        n = i-2 + (im-4)*(j-3)
        sbuf(n) = real(arr(i,j,k))
      enddo
      enddo
c
      ierr=0
      call mpi_gather(sbuf,(im-4)*(jm-4),mpi_real,
     &     rbuf,(im-4)*(jm-4),mpi_real,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
      do m = 1,np
        ipxx = mod(m-1,ipe)
        ipyy = (m-1)/ipe
!        do k = 1,km
        do j = 3,jm-2
        do i = 3,im-2
          n = i-2+(im-4)*(j-3)+(im-4)*(jm-4)*(m-1)
          ii = i-1 + lw(ipxx)
          jj = j-1 + ls(ipyy)
          if(ii.le.img .and. jj.le.jmg) then
            x2d4(ii,jj) = rbuf(n)
!            x2d4(ii,jj,k) = rbuf(n)
          endif
        enddo
        enddo
!        enddo
      enddo
        write(nunit) x2d4
      endif
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      ierr=0
      call mpi_barrier(comm,ierr)
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
      enddo
c
c
      return
      end
c
#endif

c------------------------------------------------------------------------
c
      subroutine read_r8(arr,narr,nunit)
      use param
      use mod_mpi
      implicit none

      integer nunit,narr
      real*8 arr(narr)
c
      if(ip.eq.imaster) then
        read(nunit) arr
      endif
c
      call bcast_dble(arr,narr)
c
      return
      end

