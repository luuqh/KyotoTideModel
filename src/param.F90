!
!   $Id: param.F90 10 2008-11-08 10:07:39Z ishikawa $
!
      module param

	  implicit none

#ifdef GL11M
	  integer,parameter :: imm=360,jmm=150
	  integer,parameter :: km=34,kadp=3,kshrd=7,kmlmin=1
          integer,parameter :: imn=12,nsf=360,nmfit=5,nhour=24*31
#endif
#ifdef EQ100M
      integer,parameter :: imm=10,jmm=100
      integer,parameter :: km=225,kadp=20,kshrd=200,kmlmin=20
      integer,parameter :: imn=12,nsf=360,nmfit=5,nhour=24*31
#endif
#ifdef PC68M
      integer,parameter :: imm=1070-5,jmm=609-5
      integer,parameter :: km=78,kadp=5,kshrd=20,kmlmin=3
!      integer,parameter :: imn=12,nsf=360,nmfit=5,nhour=24*31
#ifdef NCEPSF
      integer,parameter :: imn=12,nmfit=5,nhour=24*31,nsf=366
#else
      integer,parameter :: imn=12,nmfit=5,nhour=24*31,nsf=360
#endif
#endif

#ifdef JP44
      integer,parameter :: km=34,kadp=3,kshrd=7,kmlmin=1
      integer,parameter :: inc=4,jnc=4
      integer,parameter :: ig=35,jg=25
      integer,parameter :: ism=97,jsm=102
      integer,parameter :: im0=364,jm0=155
!    NOTE : indg > 0 and jndg > 0
      integer,parameter :: nrbd=711,ndpt=im0*jm0,indg=2,jndg=2
      integer,parameter :: imn=12,nsf=360,nmfit=5,nhour=24*31
#endif

#ifdef NWNPAC
      integer,parameter :: km=78,kadp=5,kshrd=20,kmlmin=3
      integer,parameter :: inc=3,jnc=3,ig=168,jg=136
!using JP68M for open bounday
      integer,parameter :: ism=135-108+4,jsm=328-288+4
      integer,parameter :: im0=265,jm0=223
!------------------------------------------
!using PC68M for open bounday
!      integer,parameter :: ism=135,jsm=328
!      integer,parameter :: im0=1070,jm0=609
!----------------------------------------------
!      integer,parameter :: nrbd=13,ndpt=im0*jm0,indg=5,jndg=5
      integer,parameter :: ndpt=im0*jm0,indg=5,jndg=5
      integer,parameter :: nrbd=6
      integer,parameter :: imn=12,nsf=366,nmfit=5,nhour=24*31
#endif

#ifdef JP68M
      integer,parameter :: km=78,kadp=5,kshrd=20,kmlmin=3
      integer,parameter :: inc=1,jnc=1,ig=258,jg=216
      integer,parameter :: ism=108,jsm=288
      integer,parameter :: im0=1070,jm0=609
      integer,parameter :: nrbd=13,ndpt=im0*jm0,indg=2,jndg=2
#ifdef NCEPSF
      integer,parameter :: imn=12,nmfit=5,nhour=24*31,nsf=366
#else
      integer,parameter :: imn=12,nmfit=5,nhour=24*31,nsf=360
#endif
#endif

#ifdef ROKKA
      integer,parameter :: km=78,kadp=5,kshrd=20,kmlmin=3
      integer,parameter :: inc=3,jnc=3,ig=126,jg=114
      integer,parameter :: ism=206,jsm=193
      integer,parameter :: im0=511,jm0=415
!*** Edit by Ishikawa-sensei 2009-02-24
!     integer,parameter :: nrbd=21,ndpt=im0*jm0,indg=5,jndg=5
!     integer,parameter :: nrbd=21,ndpt=im0*jm0,indg=20,jndg=20
      integer,parameter :: nrbd=21,ndpt=im0*jm0,indg=5,jndg=5
      integer,parameter :: imn=12,nmfit=5,nhour=24*31
#if defined(EXP2006) || defined(EXP2007)
      integer,parameter :: nsf=28*4,nws=28*8

!#elif defined(EXP2007)
!      integer,parameter :: nsf=366*8,nws=366*8
#else
      integer,parameter :: nsf=366*4
#endif
#endif
      integer :: nday

#ifdef NESTED
      integer,parameter :: imm=ig*inc+2,jmm=jg*jnc+2
#endif

#ifdef CYCLIC
      integer,parameter :: img= imm+4, jmg = jmm+5
#else
      integer,parameter :: img= imm+5, jmg = jmm+5
#endif

#if defined(ROKKA) && defined(EXP2007)
      integer,parameter :: numrm=100
#endif

#ifdef NO_MPI
      integer,parameter :: im = img, jm = jmg
      integer,parameter :: ipe = 1,jpe = 1, npe=ipe*jpe
#else

! parameters for mpi_local
#if defined(ROKKA) && defined(EXP2006) && defined(DEBUG13)
      integer,parameter :: ipe = 4,jpe = 4, npe=ipe*jpe
#else
#if defined(DEBUG20) && defined(NWNPAC)
      integer,parameter :: ipe = 2,jpe = 2, npe=ipe*jpe
#else
      integer,parameter :: ipe = 8,jpe = 4, npe=ipe*jpe
#endif
#endif

#ifdef PC68M
      integer,parameter :: ipe = 8,jpe = 8, npe=ipe*jpe
#endif

#ifdef JP44
      integer,parameter :: ipe = 4,jpe = 4, npe=ipe*jpe
#endif
      integer,parameter :: im = (imm+ipe)/ipe+4, jm = (jmm+jpe)/jpe+4
      integer,parameter :: ijmax = max0(im,jm)
#endif

! --  division of array
      integer ls(0:jpe-1),ln(0:jpe-1),ljm(0:jpe-1)
      integer lw(0:ipe-1),le(0:ipe-1),lim(0:ipe-1)
	  integer iml,jml

#ifdef BNDTIDE
      integer,parameter :: tcm = 61       ! tidal components
      integer,parameter :: ncnfile = 1000 ! max no. of restart files
#endif

      end module param



! NOTE: Explit values of parameters for our BNDTIDE
! inc = 3
! jnc = 3
! ig = 126
! jg = 114
! imm = ig*inc+2 = 380
! jmm = jg*jnc+2 = 344
! img = imm+5 = 385
! jmg = jmm+5 = 349
! ipe = 8
! jpe = 4
! npe = ipe*jpe = 32
! im = (imm+ipe)/ipe+4 = 53.1
! jm = (jmm+jpe)/jpe+4 = 92.2