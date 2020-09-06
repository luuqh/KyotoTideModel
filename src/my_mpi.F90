c
c  re-difnied mpi subroutine 
c  $Id: my_mpi.F90 12 2008-12-11 11:12:34Z ishikawa $
c-----------------------------------------------------------------------
c
      subroutine my_mpi_init()
 
      use param
      use mod_mpi
      implicit none

      integer ierr,n
c
#ifdef NO_MPI
      ip = 0
      np = 1
      imaster = 0
      ls(0) = 1
      ln(0) = jmg
      ljm(0) = jm
      jml = jm
      lw(0) = 1
      le(0) = img
      lim(0) = im
      iml = im
#else
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world, ip, ierr)
      call mpi_comm_size(mpi_comm_world, np, ierr)

      call mpi_comm_dup(mpi_comm_world, comm, ierr)

      imaster = 0

!   set 

      ls(0) = 1
      ln(0) = ls(0) + jm-1
      do n = 1,jpe-1
       ls(n) = ls(n-1)+jm-4
       ln(n) = ln(n-1)+jm-4
       ljm(n-1) = jm
      enddo

      lw(0) = 1
      le(0) = lw(0) + im-1
      do n = 1,ipe-1
       lw(n) = lw(n-1)+im-4
       le(n) = le(n-1)+im-4
       lim(n-1) = im
      enddo

!cc   no good cpu number: 40,39,37,35,33,30
!cc   good cpu number   : 38,36,34,32,31,14,10  for jmg=609
!      if(ls(np-1).ge.jmm)then
!         if(ip.eq.imaster)
!     &        write(*,*)'(my_mpi.F) illegal number of cpu!',
!     &         jmm,jm,ls(np-1)
!         stop
!      endif

      ln(jpe-1) = jmg
      ljm(jpe-1) = jm - jpe*(jm-4)+(jmg-4)
      le(ipe-1) = img
      lim(ipe-1) = im - ipe*(im-4)+(img-4)

      ip_x = mod(ip,ipe)
      ip_y = ip/ipe
      iml = lim(ip_x)
      jml = ljm(ip_y)

      if(ip_y.eq.0) then
        south_id = mpi_proc_null
      else
        south_id = ip - ipe
      endif
      if(ip_y.eq.jpe-1) then
        north_id = mpi_proc_null
      else
        north_id = ip + ipe
      endif

      if(ip_x.eq.0) then
#ifdef CYCLIC
        west_id = ipe*(ip_y+1)-1
#else
        west_id = mpi_proc_null
#endif
      else
        west_id = ip - 1
      endif
      if(ip_x.eq.ipe-1) then
#ifdef CYCLIC
        east_id = (ipe*ip_y)
#else
        east_id = mpi_proc_null
#endif
      else
        east_id = ip + 1
      endif
#ifdef DEBUG
      write(*,*) 'np =',np,' ip =',ip,' master =',imaster
     &  ,' ip_xy ',ip_x,ip_y,'  local ijm = ',iml,jml
     &  ,' south =',south_id,' north =',north_id
     &  ,' west =',west_id,' east =',east_id
     &  ,' ipe =',ipe,' jpe =',jpe
      if(ip.eq.imaster)
     &     write(*,*)'ls: ',ls,' ln: ',ln,' lw: ',lw,' le: ',le,
     &     ' ljm: ',ljm,' lim: ',lim
#endif

#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine my_mpi_exit( )
 
      use param
      use mod_mpi
      implicit none

      integer ierr
c
#ifndef NO_MPI
      call mpi_finalize(ierr)
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c Written by Luu Quang Hung on August 6th 2009
c
      subroutine bcast_char(buf,narr)

      use param
      use mod_mpi
      implicit none

      character(len=80) buf(narr)
      integer narr,ierr

#ifndef NO_MPI

      call mpi_bcast(buf,narr,mpi_character,imaster,comm,ierr)

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c Written by Luu Quang Hung on August 6th 2009
c
      subroutine bcast_char1(buf,narr)

      use param
      use mod_mpi
      implicit none

      character(len=80) buf
      integer narr,ierr

#ifndef NO_MPI

      call mpi_bcast(buf,narr,mpi_character,imaster,comm,ierr)

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine bcast_int(buf,narr)

      use param
      use mod_mpi
      implicit none

      integer buf(narr),narr,ierr

#ifndef NO_MPI

      call mpi_bcast(buf,narr,mpi_integer,imaster,comm,ierr)

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine bcast_int1(buf,narr)

      use param
      use mod_mpi
      implicit none

      integer buf,narr,ierr
c
#ifndef NO_MPI
      call mpi_bcast(buf,narr,mpi_integer,imaster,comm,ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
!-----------------------------------------------------------------------

      subroutine bcast_logical1(buf)

      use param
      use mod_mpi
      implicit none

      logical,intent(inout) :: buf
      integer narr,ierr

      narr = 1
#ifndef NO_MPI
      call mpi_bcast(buf,narr,mpi_logical,imaster,comm,ierr)

      if(ierr/=0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif

#endif
      return
      end

!-----------------------------------------------------------------------

      subroutine bcast_dble(buf,narr)

      use param
      use mod_mpi
      implicit none

      integer narr,ierr
      real*8 buf(narr)
c
#ifndef NO_MPI
      call mpi_bcast(buf,narr,mpi_double_precision,
     &   imaster,comm,ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine bcast_dble1(buf,narr)

      use param
      use mod_mpi
      implicit none

      integer narr,ierr
      real*8 buf
c
#ifndef NO_MPI
      call mpi_bcast(buf,narr,mpi_double_precision,
     &   imaster,comm,ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine bcast_real(buf,narr)

      use param
      use mod_mpi
      implicit none

      integer narr,ierr
      real buf(narr)
c
#ifndef NO_MPI
      call mpi_bcast(buf,narr,mpi_real,imaster,comm,ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine bcast_real1(buf,narr)

      use param
      use mod_mpi
      implicit none

      integer narr,ierr
      real buf
c
#ifndef NO_MPI
      call mpi_bcast(buf,narr,mpi_real,imaster,comm,ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_2d_s1(arr,ntag)

      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i
      real*8 arr(im,jm)
c
#ifndef NO_MPI
      do i = 1,im
        my_rbuf(i,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),im,mpi_double_precision,
     &    south_id,ntag,comm,req_in(ntag),ierr)
c
c
      do i = 1,im
        my_sbuf(i,ntag) = arr(i,jm-3)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),im,mpi_double_precision,
     &    north_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_2d_s1(arr,ntag)
      use param
      use mod_mpi

      implicit none
      integer ierr,ntag,i
      real*8 arr(im,jm)
c
#ifndef NO_MPI
c wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)

      do i = 1,im
        arr(i,1) = my_rbuf(i,ntag)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_2d_w1(arr,ntag)
      use param
      use mod_mpi

      implicit none
      integer ierr,ntag,j
      real*8 arr(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm
        my_rbuf(j,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),jm,mpi_double_precision,
     &    west_id,ntag,comm,req_in(ntag),ierr)
c
c
      do j = 1,jm
        my_sbuf(j,ntag) = arr(iml-3,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),jm,mpi_double_precision,
     &    east_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_2d_w1(arr,ntag)
      use param
      use mod_mpi

      implicit none
      integer ierr,ntag,j
      real*8 arr(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr(1,j) = arr(iml-3,j)
      enddo
#endif
#else
c wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)

      do j = 1,jm
        arr(1,j) = my_rbuf(j,ntag)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_2d_s1p(arr1,arr2,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do i = 1,im*2
        my_rbuf(i,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),im*2,mpi_double_precision,
     &    south_id,ntag,comm,req_in(ntag),ierr)
c
c
      do i = 1,im
        my_sbuf(i,ntag) = arr1(i,jm-3)
        my_sbuf(i+im,ntag) = arr2(i,jm-3)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),im*2,mpi_double_precision,
     &    north_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_2d_s1p(arr1,arr2,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
c wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do i = 1,im
        arr1(i,1) = my_rbuf(i,ntag)
        arr2(i,1) = my_rbuf(i+im,ntag)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_2di_s1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i
      integer arr(im,jm)
c
#ifndef NO_MPI
      do i = 1,im
        imy_rbuf(i,ntag) = 0
      enddo
c
c
      call mpi_irecv(imy_rbuf(1,ntag),im,mpi_integer,south_id,
     &    ntag,comm,req_in(ntag),ierr)
c
      do i = 1,im
        imy_sbuf(i,ntag) = arr(i,jm-3)
      enddo
c
      call mpi_isend(imy_sbuf(1,ntag),im,mpi_integer,north_id,
     &    ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_2di_s1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i
      integer arr(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do i = 1,im
        arr(i,1) = imy_rbuf(i,ntag)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_2di_w1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,j
      integer arr(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm
        imy_rbuf(j,ntag) = 0
      enddo
c
      call mpi_irecv(imy_rbuf(1,ntag),jm,mpi_integer,west_id,
     &    ntag,comm,req_in(ntag),ierr)
c
      do j = 1,jm
        imy_sbuf(j,ntag) = arr(iml-3,j)
      enddo
c
      call mpi_isend(imy_sbuf(1,ntag),jm,mpi_integer,east_id,
     &    ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_2di_w1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,j
      integer arr(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr(1,j) = arr(iml-3,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do j = 1,jm
        arr(1,j) = imy_rbuf(j,ntag)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_3d_s1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i,k,n
      real*8 arr(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)
        my_rbuf(i,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),im*(km+2),mpi_double_precision,
     &    south_id,ntag,comm,req_in(ntag),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag) = arr(i,jm-3,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),im*(km+2),mpi_double_precision,
     &    north_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_3d_s1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i,k,n
      real*8 arr(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr(i,1,k) = my_rbuf(n,ntag)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_3d_s1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,jm-3,k)
        my_sbuf(n+im*(km+2),ntag1) = arr2(i,jm-3,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_3d_s1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,1,k) = my_rbuf(n,ntag1)
        arr2(i,1,k) = my_rbuf(n+im*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_2d_n1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i
      real*8 arr(im,jm)
c
#ifndef NO_MPI
      do i = 1,im
        my_rbuf(i,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),im,mpi_double_precision,
     &    north_id,ntag,comm,req_in(ntag),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag) = arr(i,4)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),im,mpi_double_precision,
     &    south_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_2d_n1(arr,ntag)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag,i
      real*8 arr(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do i = 1,im
        arr(i,jm) = my_rbuf(i,ntag)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_2d_n1p(arr1,arr2,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do i = 1,im*2
        my_rbuf(i,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),im*2,mpi_double_precision,
     &    north_id,ntag,comm,req_in(ntag),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag) = arr1(i,4)
        my_sbuf(im+i,ntag) = arr2(i,4)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),im*2,mpi_double_precision,
     &    south_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_2d_n1p(arr1,arr2,ntag)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do i = 1,im
        arr1(i,jm) = my_rbuf(i,ntag)
        arr2(i,jm) = my_rbuf(i+im,ntag)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_3d_n1(arr,ntag)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag,i,k,n
      real*8 arr(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)
        my_rbuf(i,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),im*(km+2),mpi_double_precision,
     &    north_id,ntag,comm,req_in(ntag),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag) = arr(i,4,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),im*(km+2),mpi_double_precision,
     &    south_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_3d_n1(arr,ntag)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,i,k,n
      real*8 arr(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr(i,jm,k) = my_rbuf(n,ntag)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_n1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,3,k)
        my_sbuf(n+im*(km+2),ntag1) = arr2(i,3,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_n1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,jm-1,k) = my_rbuf(n,ntag1)
        arr2(i,jm-1,k) = my_rbuf(n+im*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_e1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(3,j,k)
        my_sbuf(n+jm*(km+2),ntag1) = arr2(3,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t3d_e1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(iml-1,j,k) = arr1(3,j,k)
        arr2(iml-1,j,k) = arr2(3,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(iml-1,j,k) = my_rbuf(n,ntag1)
        arr2(iml-1,j,k) = my_rbuf(n+jm*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_e2p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(4,j,k)
        my_sbuf(n+jm*(km+2),ntag1) = arr2(4,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t3d_e2p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(iml,j,k) = arr1(4,j,k)
        arr2(iml,j,k) = arr2(4,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(iml,j,k) = my_rbuf(n,ntag1)
        arr2(iml,j,k) = my_rbuf(n+jm*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_w1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(iml-2,j,k)
        my_sbuf(n+jm*(km+2),ntag1) = arr2(iml-2,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t3d_w1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(2,j,k) = arr1(iml-2,j,k)
        arr2(2,j,k) = arr2(iml-2,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(2,j,k) = my_rbuf(n,ntag1)
        arr2(2,j,k) = my_rbuf(n+jm*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_w2p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(iml-3,j,k)
        my_sbuf(n+jm*(km+2),ntag1) = arr2(iml-3,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2)*2,
     &  mpi_double_precision,east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t3d_w2p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,j,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(1,j,k) = arr1(iml-3,j,k)
        arr2(1,j,k) = arr2(iml-3,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(1,j,k) = my_rbuf(n,ntag1)
        arr2(1,j,k) = my_rbuf(n+jm*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_n1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,3,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_n1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,jml-1,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t3d_s1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,jml-2,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_s1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,2,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t3d_w1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(iml-2,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_w1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(2,j,k) = arr1(iml-2,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(2,j,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t3d_w2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(iml-3,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_w2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(1,j,k) = arr1(iml-3,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(1,j,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t3d_e1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(3,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_e1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(iml-1,j,k) = arr1(3,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(iml-1,j,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t3d_e2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do j = 1,jm*(km+2)
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j+jm*k
        my_sbuf(n,ntag1) = arr1(4,j,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*(km+2),
     &  mpi_double_precision,west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_e2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do k = 0,km+1
      do j = 1,jm
        arr1(iml,j,k) = arr1(4,j,k)
      enddo
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do j = 1,jm
        n = j + jm*k
        arr1(iml,j,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_n2p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,4,k)
        my_sbuf(n+im*(km+2),ntag1) = arr2(i,4,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_n2p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,jm,k) = my_rbuf(n,ntag1)
        arr2(i,jm,k) = my_rbuf(n+im*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_s1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,jm-2,k)
        my_sbuf(n+im*(km+2),ntag1) = arr2(i,jm-2,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t3d_s1p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,2,k) = my_rbuf(n,ntag1)
        arr2(i,2,k) = my_rbuf(n+im*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_s1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do i = 1,im*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*2,mpi_double_precision,
     &    south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,jm-2)
        my_sbuf(i+im,ntag1) = arr2(i,jm-2)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*2,mpi_double_precision,
     &    north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_s1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,2) = my_rbuf(i,ntag1)
        arr2(i,2) = my_rbuf(i+im,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_n1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do i = 1,im*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*2,mpi_double_precision,
     &    north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,3)
        my_sbuf(i+im,ntag1) = arr2(i,3)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*2,mpi_double_precision,
     &    south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_n1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,jm-1) = my_rbuf(i,ntag1)
        arr2(i,jm-1) = my_rbuf(i+im,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_e1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*2,mpi_double_precision,
     &    east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(3,j)
        my_sbuf(j+jm,ntag1) = arr2(3,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*2,mpi_double_precision,
     &    west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_e1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(iml-1,j) = arr1(3,j)
        arr2(iml-1,j) = arr2(3,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(iml-1,j) = my_rbuf(j,ntag1)
        arr2(iml-1,j) = my_rbuf(j+jm,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_e2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*2,mpi_double_precision,
     &    east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(4,j)
        my_sbuf(j+jm,ntag1) = arr2(4,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*2,mpi_double_precision,
     &    west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_e2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(iml,j) = arr1(4,j)
        arr2(iml,j) = arr2(4,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(iml,j) = my_rbuf(j,ntag1)
        arr2(iml,j) = my_rbuf(j+jm,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_e2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm,mpi_double_precision,
     &    east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(4,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm,mpi_double_precision,
     &    west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_e2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(iml,j) = arr1(4,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(iml,j) = my_rbuf(j,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_w1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*2,mpi_double_precision,
     &    west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(iml-2,j)
        my_sbuf(j+jm,ntag1) = arr2(iml-2,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*2,mpi_double_precision,
     &    east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_w1p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(2,j) = arr1(iml-2,j)
        arr2(2,j) = arr2(iml-2,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(2,j) = my_rbuf(j,ntag1)
        arr2(2,j) = my_rbuf(j+jm,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_w2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm*2
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm*2,mpi_double_precision,
     &    west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(iml-3,j)
        my_sbuf(j+jm,ntag1) = arr2(iml-3,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm*2,mpi_double_precision,
     &    east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_w2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(1,j) = arr1(iml-3,j)
        arr2(1,j) = arr2(iml-3,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(1,j) = my_rbuf(j,ntag1)
        arr2(1,j) = my_rbuf(j+jm,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_w2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm,mpi_double_precision,
     &    west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(iml-3,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm,mpi_double_precision,
     &    east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_w2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(1,j) = arr1(iml-3,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(1,j) = my_rbuf(j,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_w1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm,mpi_double_precision,
     &    west_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(iml-2,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm,mpi_double_precision,
     &    east_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_w1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(2,j) = arr1(iml-2,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(2,j) = my_rbuf(j,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_e1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do j = 1,jm
        my_rbuf(j,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),jm,mpi_double_precision,
     &    east_id,ntag1,comm,req_in(ntag1),ierr)
c
      do j = 1,jm
        my_sbuf(j,ntag1) = arr1(3,j)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),jm,mpi_double_precision,
     &    west_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_e1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,j
      real*8 arr1(im,jm)
c
#ifdef NO_MPI
#ifdef CYCLIC
      do j = 1,jm
        arr1(iml-1,j) = arr1(3,j)
      enddo
#endif
#else
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do j = 1,jm
        arr1(iml-1,j) = my_rbuf(j,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_s1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do i = 1,im
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im,mpi_double_precision,
     &    south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,jml-2)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im,mpi_double_precision,
     &    north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_s1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,2) = my_rbuf(i,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
cc-----------------------------------------------------------------------
c
      subroutine exch_t2d_n1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do i = 1,im
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im,mpi_double_precision,
     &    north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,3)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im,mpi_double_precision,
     &    south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_n1(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,jml-1) = my_rbuf(i,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_3d_n2(arr,ntag,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,ntag2,i,k,n
      real*8 arr(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)*2
        my_rbuf(i,ntag) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag),im*(km+2)*2,mpi_double_precision,
     &   north_id,ntag,comm,req_in(ntag),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag) = arr(i,3,k)
        my_sbuf(n+im*(km+2),ntag) = arr(i,4,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag),im*(km+2)*2,mpi_double_precision,
     &   south_id,ntag,comm,req_out(ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_3d_n2(arr,ntag,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag,ntag2,i,k,n
      real*8 arr(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag),status(1,ntag),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr(i,jm-1,k) = my_rbuf(n,ntag)
        arr(i,jm,k) = my_rbuf(n+im*(km+2),ntag)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag),status(1,ntag),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sum_all(arr,narr)
      use param
      use mod_mpi
      implicit none

      integer ierr,narr,n
      real*8 arr(narr)
      real*8 sbuf(narr),rbuf(narr)
c
#ifndef NO_MPI
      do n = 1,narr
        sbuf(n) = arr(n)
      enddo
c
      call mpi_allreduce(sbuf,rbuf,narr,mpi_double_precision,
     &   mpi_sum,comm,ierr)
c
      do n = 1,narr
        arr(n) = rbuf(n)
      enddo
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine isum_all(arr,narr)
      use param
      use mod_mpi
      implicit none

      integer ierr,narr,n
      integer arr(narr)
      integer sbuf(narr),rbuf(narr)
c
#ifndef NO_MPI
      do n = 1,narr
        sbuf(n) = arr(n)
      enddo
c
      call mpi_allreduce(sbuf,rbuf,narr,mpi_integer,
     &   mpi_sum,comm,ierr)
c
      do n = 1,narr
        arr(n) = rbuf(n)
      enddo
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sum_all1(arr,narr)
      use param
      use mod_mpi
      implicit none

      integer ierr,narr,n
      real*8 arr,sbuf,rbuf
c
#ifndef NO_MPI
        sbuf = arr
c
      call mpi_allreduce(sbuf,rbuf,narr,mpi_double_precision,
     &   mpi_sum,comm,ierr)
c
        arr = rbuf
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine isum_all1(arr,narr)
      use param
      use mod_mpi
      implicit none

      integer ierr,narr
      integer arr,sbuf,rbuf
c
#ifndef NO_MPI
        sbuf = arr
c
      call mpi_allreduce(sbuf,rbuf,narr,mpi_integer,
     &   mpi_sum,comm,ierr)
c
        arr = rbuf
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sum_mst1(arr,narr)
      use param
      use mod_mpi

      integer ierr,narr
      real*8 arr
      real*8 sbuf,rbuf
c
#ifndef NO_MPI
        sbuf = arr

      call mpi_reduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_sum,imaster,comm,ierr)

      if(ip.eq.imaster) then
          arr = rbuf
      endif

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sum_mst(arr,narr)
      use param
      use mod_mpi
      implicit none

      integer ierr,narr,n
      real*8 arr(narr)
      real*8 sbuf(narr),rbuf(narr)
c
#ifndef NO_MPI
      do n = 1,narr
        sbuf(n) = arr(n)
      enddo
c
      call mpi_reduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_sum,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
        do n = 1,narr
          arr(n) = rbuf(n)
        enddo
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
c-----------------------------------------------------------------------
c
      subroutine max_mst1(arr,narr)
      use param
      use mod_mpi

      integer ierr,narr
      real*8 arr,sbuf,rbuf
c
#ifndef NO_MPI
        sbuf = arr

      call mpi_reduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_max,imaster,comm,ierr)

      if(ip.eq.imaster) then
          arr = rbuf
      endif

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine max_mst(arr,narr)

      use param
      use mod_mpi
      implicit none

      integer ierr,narr,n
      real*8 arr(narr),sbuf(narr),rbuf(narr)
c
#ifndef NO_MPI
      do n = 1,narr
        sbuf(n) = arr(n)
      enddo
c
      call mpi_reduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_max,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
        do n = 1,narr
          arr(n) = rbuf(n)
        enddo
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
c-----------------------------------------------------------------------
c
      subroutine max_all(arr,narr)
      use param
      use mod_mpi
      implicit none
      integer ierr,narr,n
      real*8 arr(narr),sbuf(narr),rbuf(narr)
c
#ifndef NO_MPI
      do n = 1,narr
        sbuf(n) = arr(n)
      enddo
c
      call mpi_allreduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_max,comm,ierr)

        do n = 1,narr
          arr(n) = rbuf(n)
        enddo

      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine min_all(arr,narr)
      use param
      use mod_mpi
      implicit none
      integer ierr,narr,n
      real*8 arr(narr),sbuf(narr),rbuf(narr)
c
#ifndef NO_MPI
      do n = 1,narr
        sbuf(n) = arr(n)
      enddo
c
      call mpi_allreduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_min,comm,ierr)
c
        do n = 1,narr
          arr(n) = rbuf(n)
        enddo
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine max_all1(arr,narr)
      use param
      use mod_mpi
      implicit none
      integer ierr,narr,n
      real*8 arr,sbuf,rbuf
c
#ifndef NO_MPI
        sbuf = arr

      call mpi_allreduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_max,comm,ierr)

          arr = rbuf
          
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine min_all1(arr,narr)
      use param
      use mod_mpi
      implicit none

      integer ierr,narr,n
      real*8 arr,sbuf,rbuf
c
#ifndef NO_MPI
        sbuf = arr
c
      call mpi_allreduce(sbuf,rbuf,narr,mpi_double_precision,
     &    mpi_min,comm,ierr)
c
          arr = rbuf
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine imax_mst(arr,narr)
      use param
      use mod_mpi
      implicit none

      integer ierr,narr,arr(narr),sbuf(narr),rbuf(narr),n
c
#ifndef NO_MPI
      do n = 1,narr
        sbuf(n) = arr(n)
      enddo
c
      call mpi_reduce(sbuf,rbuf,narr,mpi_integer,mpi_max,imaster
     &    ,comm,ierr)
c
      if(ip.eq.imaster) then
        do n = 1,narr
          arr(n) = rbuf(n)
        enddo
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
c------------------------------------------------------------------------
c
      subroutine max_at_2d(amax,iamax,jamax,nm)
      use param
      use mod_mpi
      implicit none

      integer ierr,iamax(nm),jamax(nm),itag,n,nm
      real*8 amax(nm),sbuf(2,nm),rbuf(2,nm)
c
#ifndef NO_MPI
      do n = 1,nm
        sbuf(1,n) = amax(n)
        sbuf(2,n) = dble(img*(ls(ip_y)+jamax(n)-1)+
     &     lw(ip_x)-1+iamax(n))
      enddo
c
      call mpi_reduce(sbuf,rbuf,nm,mpi_2double_precision,
     &   mpi_maxloc,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
        do n = 1,nm
          amax(n) = rbuf(1,n)
          jamax(n) = nint(rbuf(2,n))/img
          iamax(n) = nint(rbuf(2,n)) - img*(jamax(n))
        enddo
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
c
c------------------------------------------------------------------------
c
      subroutine min_at_2d(amin,iamin,jamin,nm)
      use param
      use mod_mpi
      implicit none

      integer ierr,iamin(nm),jamin(nm),itag,n,nm
      real*8 amin(nm),sbuf(2,nm),rbuf(2,nm)
c
#ifndef NO_MPI
      do n = 1,nm
        sbuf(1,n) = amin(n)
        sbuf(2,n) = dble(img*(ls(ip_y)+jamin(n)-1)
     &     +lw(ip_x)-1+iamin(n))
      enddo
c
      call mpi_reduce(sbuf,rbuf,nm,mpi_2double_precision,
     &   mpi_minloc,imaster,comm,ierr)
c
      if(ip.eq.imaster) then
        do n = 1,nm
          amin(n) = rbuf(1,n)
          jamin(n) = nint(rbuf(2,n))/img
          iamin(n) = nint(rbuf(2,n)) - img*(jamin(n))
        enddo
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
c-----------------------------------------------------------------------
c
      subroutine exch_t3d_n2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,4,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_n2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,jm,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
cc-----------------------------------------------------------------------
c
      subroutine exch_t3d_s2p(arr1,arr2,ntag1,ntag2)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,jm-3,k)
        my_sbuf(n+im*(km+2),ntag1) = arr2(i,jm-3,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2)*2,
     &  mpi_double_precision,north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t3d_s2p(arr1,arr2,ntag1,ntag2)

      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,ntag2,i,k,n
      real*8 arr1(im,jm,0:km+1),arr2(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,1,k) = my_rbuf(n,ntag1)
        arr2(i,1,k) = my_rbuf(n+im*(km+2),ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_s2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do i = 1,im*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*2,mpi_double_precision,
     &    south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,jm-3)
        my_sbuf(i+im,ntag1) = arr2(i,jm-3)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*2,mpi_double_precision,
     &    north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_s2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,1) = my_rbuf(i,ntag1)
        arr2(i,1) = my_rbuf(i+im,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_s2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none
      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do i = 1,im
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im,mpi_double_precision,
     &    south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,jml-3)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im,mpi_double_precision,
     &    north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_s2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,1) = my_rbuf(i,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
cc-----------------------------------------------------------------------
c
      subroutine exch_t3d_s2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
      do i = 1,im*(km+2)
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,south_id,ntag1,comm,req_in(ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i+im*k
        my_sbuf(n,ntag1) = arr1(i,jm-3,k)
      enddo
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*(km+2),
     &  mpi_double_precision,north_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c

c-----------------------------------------------------------------------
c
      subroutine wait_t3d_s2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i,k,n
      real*8 arr1(im,jm,0:km+1)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do k = 0,km+1
      do i = 1,im
        n = i + im*k
        arr1(i,1,k) = my_rbuf(n,ntag1)
      enddo
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_n2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
      do i = 1,im*2
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im*2,mpi_double_precision,
     &    north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,4)
        my_sbuf(i+im,ntag1) = arr2(i,4)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im*2,mpi_double_precision,
     &    south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_n2p(arr1,arr2,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm),arr2(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,jm) = my_rbuf(i,ntag1)
        arr2(i,jm) = my_rbuf(i+im,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine exch_t2d_n2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
      do i = 1,im
        my_rbuf(i,ntag1) = 0.
      enddo
c
      call mpi_irecv(my_rbuf(1,ntag1),im,mpi_double_precision,
     &    north_id,ntag1,comm,req_in(ntag1),ierr)
c
      do i = 1,im
        my_sbuf(i,ntag1) = arr1(i,4)
      enddo
c
      call mpi_isend(my_sbuf(1,ntag1),im,mpi_double_precision,
     &    south_id,ntag1,comm,req_out(ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine wait_t2d_n2(arr1,ntag1)
      use param
      use mod_mpi
      implicit none

      integer ierr,ntag1,i
      real*8 arr1(im,jm)
c
#ifndef NO_MPI
c  wait for receive
      call mpi_wait(req_in(ntag1),status(1,ntag1),ierr)
c
      do i = 1,im
        arr1(i,jml) = my_rbuf(i,ntag1)
      enddo
c
c  wait for send
      call mpi_wait(req_out(ntag1),status(1,ntag1),ierr)
c
      if(ierr.ne.0) then
        write(*,*) 'ip=',ip,'  ierr =',ierr
      endif
c
#endif
      return
      end
c
c-----------------------------------------------------------------------
