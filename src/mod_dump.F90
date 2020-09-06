!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   IMPORT RESULT FROM NAOVARJ2 TO ASSIMILATION SYSTEM
!
!   4DVAR output code :  Luu Quang Hung
!   Update            :  2010.07.02
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
! $ main - Luu Quang Hung
!-------------------------------------------------------------------

      MODULE mod_dump

      use param
      implicit none

      private
      public caldump

      contains

! $ caldump - Luu Quang Hung
!-------------------------------------------------------------------
      SUBROUTINE caldump
!-------------------------------------------------------------------

      use mod_mpi
      implicit none

#include "common.h"

      real*8 dg(img,jmg)
      real*8 kex,key,alpha,beta,gama,fdump
      real*8 f0,f1,f2,a,b
      integer*4 i0,j0,i1,j1,i,j,dumpcase
	  
      dumpcase=6
	  
      if(ip==imaster) then
	    
         select case (dumpcase)

         case (0)
		 
            do j = 1,jmg
               do i = 1,img
                 dg(i,j) = 5.5d-6
               enddo
            enddo	
       	  
         case (1) ! constant Rok + exponetial dumping SW
		 
            fdump = .5d-6
            i0 = 24
            j0 = 1
            i1 = 40
            j1 = 50
            alpha = 30.  ! peak 10 times larger than average
            beta = 1./10  ! at i1, dump is 1/10. larger than at origin i0 
            gama = 1./10. ! at j1, dump is 1/10. larger than at origin j0
            ! note that: alpha*beta>1, alpha*gama>1
            do j = 1,jmg
               do i = 1,img
                  kex=log((beta*alpha-1)/(alpha-1))*((i-i0)/(i1-i0))**2
                  key=log((gama*alpha-1)/(alpha-1))*((j-j0)/(j1-j0))**2
                  dg(i,j) = fdump*(1+(alpha-1)*exp(kex+key))
               enddo
            enddo
		 
         case (2) ! constant Rok + linear dumping SW
		 
            f0 = 4.0d-6
            f1 = 5.5d-6 
            i0 = 110
            j0 = 150
            j1 = 200
            do j = 1,jmg
               do i = 1,img
                  dg(i,j) = f0
               enddo
            enddo
            a = (f0-f1)/(j1-j0)
            b = f0-a*j1
            do i = 1,i0			
               do j = 1,j0                    			   
                  dg(i,j) = f1
               enddo
               do j = j0+1,j1
                  dg(i,j) = a*j+b
               enddo
            enddo			

         case (3)  ! constant + sponze zone Rok
		 
            f0 = 20.0d-6
            f1 = 2.d-6
            i1 = 50 
            j1 = 50
            do j = 1,jmg
               do i = 1,img
                  dg(i,j) = f0
               enddo
            enddo
            do j = 1+j1,jmg-j1
               do i = 1+i1,img-j1
                  dg(i,j) = f1
               enddo
            enddo

         case (4) ! constant + sponze zone Rok + linear dumping SW      
		 
            f0 = 20.0d-6
            f1 = 1.d-6
            i1 = 50 
            j1 = 50
            do j = 1,jmg
               do i = 1,img
                  dg(i,j) = f0
               enddo
            enddo
            do j = 1+j1,jmg-j1
               do i = 1+i1,img-j1
                  dg(i,j) = f1
               enddo
            enddo		 
            i0 = 140
            j0 = 100
            do i = 1,i0			
               do j = 1,j0                    			   
                  dg(i,j) = f0
               enddo
            enddo	
		 
         case (5)  ! constant + SW dumping
		 
            f0 = 5.d-3
            f1 = 5.d-5
            f2 = 5.d-7
            i1 = 20 
            j1 = 20
            do j = 1,jmg
               do i = 1,img
                  dg(i,j) = f0
               enddo
            enddo
            do j = 1+j1,jmg-j1
               do i = 1+i1,img-i1
                  dg(i,j) = f1
               enddo
            enddo		 
            do j = 1+2*j1,jmg-2*j1
               do i = 1+2*i1,img-2*i1
                  dg(i,j) = f2
               enddo
            enddo		 
            i0 = 140
            j0 = 100
            do i = 1,i0			
               do j = 1,j0                    			   
                  dg(i,j) = f0
               enddo
            enddo	

		 
         case (6) ! present
		 
            f0 = 5.d-2
            f1 = 5.d-4
            f2 = 5.d-8
            i1 = 10 
            j1 = 10
            do j = 1,jmg
               do i = 1,img
                  dg(i,j) = f0
               enddo
            enddo
            do j = 1+j1,jmg-j1
               do i = 1+i1,img-i1
                  dg(i,j) = f1
               enddo
            enddo		 
            do j = 1+2*j1,jmg-2*j1
               do i = 1+2*i1,img-2*i1
                  dg(i,j) = f2
               enddo
            enddo	 
			
         case default   ! default case
		 
            do j = 1,jmg
               do i = 1,img
                 dg(i,j) = 5.5d-6
               enddo
            enddo	
		 
         end select
		 
      endif !imaster

      call bcast_dble(dg,img*jmg)
		
      do j = 1,jml
         do i = 1,iml
            dump1(i,j) = dg(lw(ip_x)-1+i,ls(ip_y)-1+j)
            dump2(i,j) = 1.d-5
            dump3(i,j) = 0.
         enddo
      enddo			


      return
      END SUBROUTINE caldump



      END MODULE mod_dump
c ----------------------< End of program >----------------------
c \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/


