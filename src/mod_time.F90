!
!  parameters and variables for time control
!  $Id: mod_time.F90 10 2008-11-08 10:07:39Z ishikawa $
!----------------------------------------------------------------------
      module mod_time

      use datetime
      implicit none

      integer step_3d,ns_year,ns_month,ns_day,ns_hour,ns_min,ns_sec ! to make step_3d
      integer nxmin,nxhour,nxday,nxmonth,nxyear ! time steps of hour,day,month,year
      real nxsec !add by LQH

      type(TDate),save :: date_now,date_start
      type(Ttimediff),save :: diff_dtuv,diff_dttr

! jobp parameters
      logical nfirst,set_start


      end module mod_time