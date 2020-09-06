!
!  parameters for mpi
!  $Id: mod_mpi.F90 2 2007-11-15 13:07:08Z ishikawa $
!----------------------------------------------------------------------
      module mod_mpi

	  use param
	  
	  implicit none
#ifndef NO_MPI
      include (mpif.h)
#endif

!
!-- define mpp commons
      integer ip, np, comm,ip_x,ip_y,imaster
! rank,size,commnunicator
      integer south_id, north_id,east_id,west_id

#ifndef NO_MPI
!-- communication  request
      integer,parameter :: num_msg = 120
      integer status(MPI_STATUS_SIZE,num_msg)
      integer req_in(num_msg),req_out(num_msg)
	  integer imy_sbuf(ijmax*(km+2),num_msg)
	  integer imy_rbuf(ijmax*(km+2),num_msg)
      real*8 my_sbuf(ijmax*(km+2)*2,num_msg)
	  real*8 my_rbuf(ijmax*(km+2)*2,num_msg)
#endif
	  end module mod_mpi