!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                   c
!     subroutine calpp                                              c
!
!   $Id: calpp.F90 2 2007-11-15 13:07:08Z ishikawa $
!                                                                   c
!     Calculation of horizontal-mean pressure and coefficients in   c
!     the equation of state                                         c
!                                                                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calpp

      use param
!      use mod_mpi
      implicit none

#include "common.h"

      integer i,j,k

      do k = 0,km+1
        dmn(k) = 0.d0
      enddo

      do k = 1,km
      do j = 3,jml-2
      do i = 3,iml-2
        dmn(k)=dmn(k)+rho(i,j,k)*volt(i,j,k)
      enddo
      enddo
      enddo

      call sum_all(dmn,km+2)

      ddmna = 0.d0

      do k=1,km
        ddmna=ddmna+dmn(k)/ttvol
        dmn(k)=dmn(k)/tvol(k)
      enddo
      pd(1)=(9.81d-4)*dzz(1)*(dmn(1)+1.d0)
      pm(1)=0.d0
      do k=2,km
        pd(k)=pd(k-1)+(4.905d-4)*dzz(k)*(dmn(k-1)+dmn(k)+2.d0)
        pm(k)=pm(k-1)+(9.81d-4)*dz(k-1)*(dmn(k-1)+1.d0)
      enddo
      pm(km+1)=pm(km)+(9.81d-4)*dz(km)*(dmn(km)+1.d0)
      ddmna=ddmna+1.d0
      ddmnar=1.d0/ddmna

      return
      end
