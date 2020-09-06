! date managing module

!#define DEFAULT_REAL_BYTE 4
#define DEFAULT_REAL_BYTE 8

module datetime
 
  implicit none

  private

  public TDate, TTimediff, TClimdate,&
       init_date1, init_date2, init_date3,&
       init_climdate1, init_climdate2, init_climdate3,&
       timediff, year_day, year_sec, abs_timediff,replace_date,&
       replace_climdate,base_year, days_a_year_clim,&
       operator(-),operator(+),operator(*),operator(/),assignment(=)
#ifndef NO_MPI
  public bcast_Tdate

#endif

!  days_a_year_clim = 360, 365 or 366

  integer,parameter :: base_year=1950,days_a_year_clim=365

  integer,parameter :: hour_a_day=24, min_an_hour=60, sec_a_min=60,&
       min_a_day=hour_a_day*min_an_hour,&
       sec_a_day=min_a_day*sec_a_min,&
       sec_an_hour=sec_a_min*min_an_hour,days_a_week=7


  type TDate
     real(8) :: tot_sec=0.d0
     integer :: year=base_year,month=1,day=1,hour=0,min=0
     real(DEFAULT_REAL_BYTE) :: sec=0.d0
  end type TDate

  type TTimediff
     real(8) :: diff_sec=0.d0
     integer :: min=0,hour=0,day=0,week=0
     real(DEFAULT_REAL_BYTE) :: sec=0.d0
  end type TTimediff

  type TClimdate
     real(8) :: tot_sec=0.d0
     integer :: month=1,day=1,hour=0,min=0,sum_year=0
     real(DEFAULT_REAL_BYTE) :: sec=0.d0
  end type TClimdate

  interface operator(-)
     module procedure date_minus_date
     module procedure date_minus_timediff
     module procedure climdate_minus_climdate
     module procedure climdate_minus_timediff
     module procedure timediff_minus_timediff
     module procedure minus_timediff
  end interface

  interface operator(+)
     module procedure date_plus_timediff
     module procedure timediff_plus_date
     module procedure climdate_plus_climdate
     module procedure climdate_plus_timediff
     module procedure timediff_plus_climdate
     module procedure timediff_plus_timediff
  end interface

  interface operator(*)
     module procedure timediff_times_single
     module procedure timediff_times_double
     module procedure timediff_times_int
     module procedure single_times_timediff
     module procedure double_times_timediff
     module procedure int_times_timediff
     module procedure climdate_times_double
     module procedure climdate_times_single
     module procedure climdate_times_int
     module procedure double_times_climdate
     module procedure single_times_climdate
     module procedure int_times_climdate
  end interface

  interface operator(/)
     module procedure timediff_div_timediff
     module procedure timediff_div_single
     module procedure timediff_div_double
     module procedure timediff_div_int
     module procedure climdate_div_double
     module procedure climdate_div_single
     module procedure climdate_div_int
  end interface

  interface assignment(=)
     module procedure init_date4
     module procedure init_climdate4
     module procedure assign_timediff
  end interface

  interface init_date3
     module procedure init_date3_double
     module procedure init_date3_int
  end interface

  interface init_climdate3
     module procedure init_climdate3_double
     module procedure init_climdate3_int
  end interface

contains

!---------------------- initializing functions --------------------------

!------------ for date
  function init_date1(year,month,day,hour,min,sec) result(date)

    implicit none
    type(TDate) :: date
    integer,intent(in) :: year
    integer,intent(in),optional :: month,day,hour,min
    real(DEFAULT_REAL_BYTE),intent(in),optional :: sec
#ifdef DEBUG
    integer :: mday(12)
#endif

#ifdef DEBUG
    if(year<base_year) then
       write(0,*) 'error: year < base_year'
       stop
    endif
#endif
    date%year = year

    if (present(month)) then
#ifdef DEBUG
    if(month<1 .or. month>12) then
       write(0,*) 'error: month is out of range'
       stop
    endif
#endif
       date%month = month
    endif

    if (present(day)) then
#ifdef DEBUG
       mday=days_month(year)
       if(day<0 .or. day>mday(date%month)) then
          write(0,*) 'error: day is out of range'
          stop
       endif
#endif
       date%day = day
    endif

    if (present(hour)) then
#ifdef DEBUG
       if(hour<0 .or. hour>23) then
          write(0,*) 'error: hour is out of range'
          stop
       endif
#endif
       date%hour = hour
    endif
    
    if (present(min)) then
#ifdef DEBUG
       if(min<0 .or. min>59) then
          write(0,*) 'error: min is out of range'
          stop
       endif
#endif
       date%min = min
    endif

    if (present(sec)) then
#ifdef DEBUG
       if(sec<0 .or. sec>59) then
          write(0,*) 'error: sec is out of range'
          stop
       endif
#endif
       date%sec = sec
    endif

    date%tot_sec = cal_tot_sec(date)

  end function init_date1


  function init_date2(year,year_day,hour,min,sec) result(date)

    implicit none
    type(TDate) :: date
    integer,intent(in) :: year,year_day
    integer,intent(in),optional :: hour,min
    real(DEFAULT_REAL_BYTE),intent(in),optional :: sec
    
#ifdef DEBUG
    if(year<base_year) then
       write(0,*) 'error: year < base_year'
       stop
    endif
#endif

    date = year_day_to_date(year,year_day)
    if (present(hour)) then
#ifdef DEBUG
       if(hour<0 .or. hour>23) then
          write(0,*) 'error: hour is out of range'
          stop
       endif
#endif
       date%hour = hour
    endif

    if (present(min)) then
#ifdef DEBUG
       if(min<0 .or. min>59) then
          write(0,*) 'error: min is out of range'
          stop
       endif
#endif
       date%min = min
    endif

    if (present(sec)) then
#ifdef DEBUG
       if(sec<0 .or. sec>=60) then
          write(0,*) 'error: sec is out of range'
          stop
       endif
#endif
       date%sec = sec
    endif

    date%tot_sec = cal_tot_sec(date)

  end function init_date2


  function init_date3_double(year,year_sec) result(date)

    implicit none
    type(TDate) :: date
    integer,intent(in) :: year
    real(8),intent(in) :: year_sec

#ifdef DEBUG
    if(year < base_year) then
       write(0,*) 'error: year < base_year'
       stop
    endif
#endif

    date = year_sec_to_date(year,year_sec)

  end function init_date3_double

  function init_date3_int(year,year_sec) result(date)

    implicit none
    type(TDate) :: date
    integer,intent(in) :: year
    integer,intent(in) :: year_sec

#ifdef DEBUG
    if(year < base_year) then
       write(0,*) 'error: year < base_year'
       stop
    endif
    if (year_sec < 0) then
       write(0,*) 'error: year_sec < 0'
       stop
    endif
#endif

    date = year_sec_to_date(year,dble(year_sec))

  end function init_date3_int


  subroutine init_date4(date,tot_sec)

    implicit none
    type(TDate),intent(out) :: date
    real(8),intent(in) :: tot_sec

#ifdef DEBUG
    if (tot_sec < 0) then
       write(0,*) 'error: tot_sec < 0'
       stop
    endif
#endif

    date%tot_sec = tot_sec

    call update_date(date)

  end subroutine init_date4


  subroutine update_date(date)
    implicit none
    type(TDate),intent(inout) :: date
    integer :: tot_day,sum_day(2)
    integer :: cyear
    real(8) :: ysec

    tot_day = int(date%tot_sec/dble(sec_a_day))

#ifdef DEBUG
    if (date%tot_sec < 0) then
       write(0,*) 'error: tot_sec < 0'
       stop
    endif
#endif

    cyear = base_year
    sum_day(1) = 0
    sum_day(2) = days_year(cyear)
    do 
       if (tot_day >= sum_day(1) .and. tot_day < sum_day(2)) exit
       sum_day(1) = sum_day(2)
       cyear = cyear + 1
       sum_day(2) = sum_day(2) + days_year(cyear)
    end do

    ysec = date%tot_sec-dble(sum_day(1))*dble(sec_a_day)

    date = year_sec_to_date(cyear,ysec)

  end subroutine update_date

! ---------- for climdate
  function init_climdate1(month,day,hour,min,sec,sum_year) result(date)

    implicit none
    type(TClimdate) :: date
    integer,intent(in),optional :: sum_year,month,day,hour,min
    real(DEFAULT_REAL_BYTE),intent(in),optional :: sec
#ifdef DEBUG
    integer :: mday(12)
#endif

    if(present(sum_year)) date%sum_year = sum_year

    if (present(month)) then
#ifdef DEBUG
    if(month<1 .or. month>12) then
       write(0,*) 'error: month is out of range'
       stop
    endif
#endif
       date%month = month
    endif

    if (present(day)) then
#ifdef DEBUG
       mday=days_month_clim()
       if(day<0 .or. day>mday(date%month)) then
          write(0,*) 'error: day is out of range'
          stop
       endif
#endif
       date%day = day
    endif

    if (present(hour)) then
#ifdef DEBUG
       if(hour<0 .or. hour>23) then
          write(0,*) 'error: hour is out of range'
          stop
       endif
#endif
       date%hour = hour
    endif
    
    if (present(min)) then
#ifdef DEBUG
       if(min<0 .or. min>59) then
          write(0,*) 'error: min is out of range'
          stop
       endif
#endif
       date%min = min
    endif

    if (present(sec)) then
#ifdef DEBUG
       if(sec<0 .or. sec>59) then
          write(0,*) 'error: sec is out of range'
          stop
       endif
#endif
       date%sec = sec
    endif

    date%tot_sec = cal_tot_sec_clim(date)

  end function init_climdate1


  function init_climdate2(year_day,hour,min,sec,sum_year) result(date)

    implicit none
    type(TClimdate) :: date
    integer,intent(in) :: year_day
    integer,intent(in),optional :: hour,min,sum_year
    real(DEFAULT_REAL_BYTE),intent(in),optional :: sec


    if(present(sum_year)) then
       date = year_day_to_climdate(year_day,sum_year)
    else
       date = year_day_to_climdate(year_day,0)
    endif

    if (present(hour)) then
#ifdef DEBUG
       if(hour<0 .or. hour>23) then
          write(0,*) 'error: hour is out of range'
          stop
       endif
#endif
       date%hour = hour
    endif

    if (present(min)) then
#ifdef DEBUG
       if(min<0 .or. min>59) then
          write(0,*) 'error: min is out of range'
          stop
       endif
#endif
       date%min = min
    endif

    if (present(sec)) then
#ifdef DEBUG
       if(sec<0 .or. sec>=60) then
          write(0,*) 'error: sec is out of range'
          stop
       endif
#endif
       date%sec = sec
    endif

    date%tot_sec = cal_tot_sec_clim(date)

  end function init_climdate2

  function init_climdate3_double(year_sec,sum_year) result(date)

    implicit none
    type(TClimdate) :: date
    integer,intent(in),optional :: sum_year
    real(8),intent(in) :: year_sec

    if (present(sum_year)) then
       date = year_sec_to_climdate(year_sec,sum_year)
    else
       date = year_sec_to_climdate(year_sec,0)
    endif

  end function init_climdate3_double

  function init_climdate3_int(year_sec,sum_year) result(date)

    implicit none
    type(TClimdate) :: date
    integer,intent(in),optional :: sum_year
    integer,intent(in) :: year_sec

    if (present(sum_year)) then
       date = year_sec_to_climdate(dble(year_sec),sum_year)
    else
       date = year_sec_to_climdate(dble(year_sec),0)
    endif

  end function init_climdate3_int

  subroutine init_climdate4(date,tot_sec)
    implicit none
    type(TClimdate),intent(out) :: date
    real(8),intent(in) :: tot_sec

    date%tot_sec = tot_sec

    call update_climdate(date)

  end subroutine init_climdate4


  subroutine update_climdate(date)
    implicit none
    type(TClimdate),intent(inout) :: date
    integer :: tot_day,sum_day(2)
    real(8) :: ysec


    tot_day = int(date%tot_sec/dble(sec_a_day))

#ifdef DEBUG
    if (date%tot_sec < 0) then
       write(0,*) 'error: tot_sec < 0'
       stop
    endif
#endif

    date%sum_year = int(date%tot_sec/dble(days_a_year_clim)/dble(sec_a_day))
    ysec = date%tot_sec - dble(date%sum_year)*dble(days_a_year_clim)&
         *dble(sec_a_day)

    date = year_sec_to_climdate(ysec,date%sum_year)

  end subroutine update_climdate

! ---------- for timediff

  function timediff(diff_sec,week,day,hour,min,sec) 
    implicit none
    integer,intent(in),optional :: week,day,hour,min
    real(DEFAULT_REAL_BYTE),intent(in),optional ::sec
    real(8),intent(in),optional :: diff_sec
    type(TTimediff) :: timediff

    if (present(diff_sec)) then
       timediff%diff_sec = diff_sec
#ifdef DEBUG
       if(present(week)) write(0,*) 'waring: dummy argument "week" is ignored'
       if(present(day))  write(0,*) 'waring: dummy argument "day" is ignored'
       if(present(hour))  write(0,*) 'waring: dummy argument "hour" is ignored'
       if(present(min))  write(0,*) 'waring: dummy argument "min" is ignored'
       if(present(sec))  write(0,*) 'waring: dummy argument "sec" is ignored'
#endif
    else
       if(present(week)) timediff%week = week
       if(present(day)) timediff%day = day
       if(present(hour)) timediff%hour = hour
       if(present(min)) timediff%min = min
       if(present(sec)) timediff%sec = sec

       timediff%diff_sec = dble(timediff%week*days_a_week*sec_a_day)&
            +dble(timediff%day*sec_a_day)&
            +dble(timediff%hour*sec_an_hour)&
            +dble(timediff%min*sec_a_min)&
            +dble(timediff%sec)
    end if

    call update_timediff(timediff)

  end function timediff

  subroutine assign_timediff(timediff,diff_sec)
    implicit none
    type(TTimediff),intent(out) :: timediff
    real(8),intent(in) :: diff_sec


    timediff%diff_sec = diff_sec

    call update_timediff(timediff)

  end subroutine assign_timediff


  subroutine update_timediff(timediff)
    implicit none
    type(TTimediff),intent(inout) :: timediff

    call wdhms_from_sec(timediff%diff_sec,timediff%week,timediff%day,&
         timediff%hour,timediff%min,timediff%sec)

  end subroutine update_timediff

!------------------------- operating functions -------------------

! -------------  for TDate 

!  new_date = date + timediff
  function date_plus_timediff(date,timediff) result(new_date)
    implicit none
    type(TDate),intent(in) :: date
    type(TTimediff),intent(in) :: timediff
    type(TDate) :: new_date

    new_date%tot_sec = date%tot_sec + timediff%diff_sec 

    call update_date(new_date)

  end function date_plus_timediff

  function timediff_plus_date(timediff,date) result(new_date)
    implicit none
    type(TDate),intent(in) :: date
    type(TTimediff),intent(in) :: timediff
    type(TDate) :: new_date

    new_date%tot_sec = date%tot_sec + timediff%diff_sec 

    call update_date(new_date)

  end function timediff_plus_date


! new_timediff = timediff1 + timediff2
  function timediff_plus_timediff(timediff1,timediff2) result(new_timediff)
    type(TTimediff),intent(in) :: timediff1,timediff2
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff1%diff_sec + timediff2%diff_sec

    call update_timediff(new_timediff)

  end function timediff_plus_timediff

! new_timediff = timediff1 - timediff2
  function timediff_minus_timediff(timediff1,timediff2) result(new_timediff)
    type(TTimediff),intent(in) :: timediff1,timediff2
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff1%diff_sec - timediff2%diff_sec

    call update_timediff(new_timediff)

  end function timediff_minus_timediff

! new_date = date - timediff
  function date_minus_timediff(date,timediff) result(new_date)
    implicit none
    type(TDate),intent(in) :: date
    type(TTimediff),intent(in) :: timediff
    type(TDate) :: new_date

    new_date%tot_sec = date%tot_sec - timediff%diff_sec 

#ifdef DEBUG
    if (new_date%tot_sec < 0) then
       write(0,*) 'error: tot_sec < 0'
       stop
    endif
#endif

    call update_date(new_date)

  end function date_minus_timediff


!  timediff = date1 - date2
  function date_minus_date(date1,date2) result(timediff)
    implicit none
    type(TDate),intent(in) :: date1,date2
    type(TTimediff) :: timediff

    timediff%diff_sec = date1%tot_sec - date2%tot_sec

    call update_timediff(timediff)

  end function date_minus_date
  

! quot = timediff1 / timediff2
  function timediff_div_timediff(timediff1,timediff2) result(quot)
    implicit none
    type(TTimediff),intent(in) :: timediff1,timediff2
    real(DEFAULT_REAL_BYTE) :: quot

    quot = timediff1%diff_sec / timediff2%diff_sec

  end function timediff_div_timediff

! new_timediff = timediff * val
  function timediff_times_double(timediff,val) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    real(8),intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec * val

    call update_timediff(new_timediff)

  end function timediff_times_double

  function timediff_times_single(timediff,val) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    real,intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec * dble(val)

    call update_timediff(new_timediff)

  end function timediff_times_single
  
  function timediff_times_int(timediff,val) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    integer,intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec * dble(val)

    call update_timediff(new_timediff)
    
  end function timediff_times_int

! new_timediff = val * timediff 
  function double_times_timediff(val,timediff) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    real(8),intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec * val

    call update_timediff(new_timediff)

  end function double_times_timediff

  function single_times_timediff(val,timediff) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    real,intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec * dble(val)

    call update_timediff(new_timediff)

  end function single_times_timediff
  
  function int_times_timediff(val,timediff) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    integer,intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec * dble(val)

    call update_timediff(new_timediff)
    
  end function int_times_timediff

! new_timediff = timediff / val
  function timediff_div_double(timediff,val) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    real(8),intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec / val
    
    call update_timediff(new_timediff)

  end function timediff_div_double

  function timediff_div_single(timediff,val) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    real,intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec / dble(val)
    
    call update_timediff(new_timediff)

  end function timediff_div_single

  function timediff_div_int(timediff,val) result(new_timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    integer,intent(in) :: val
    type(TTimediff) :: new_timediff

    new_timediff%diff_sec = timediff%diff_sec / dble(val)
    
    call update_timediff(new_timediff)

  end function timediff_div_int


!  minus_timediff = -timediff
  function minus_timediff(timediff)
    implicit none
    type(TTimediff),intent(in) :: timediff
    type(TTimediff) :: minus_timediff

    minus_timediff%diff_sec = -timediff%diff_sec
    minus_timediff%week = -timediff%week
    minus_timediff%day = -timediff%day
    minus_timediff%hour = -timediff%hour
    minus_timediff%min = -timediff%min
    minus_timediff%sec = -timediff%sec

  end function minus_timediff

! absolute value of timediff
  function abs_timediff(timediff) 
    implicit none
    type(TTimediff),intent(in) :: timediff
    type(TTimediff) :: abs_timediff

    abs_timediff%diff_sec = abs(timediff%diff_sec)

    call update_timediff(abs_timediff)

  end function abs_timediff

! replace date member

  function replace_date(date,year,month,day,hour,min,sec) result(new_date)
    implicit none
    integer,intent(in),optional :: year,month,day,hour,min
    real(DEFAULT_REAL_BYTE),intent(in),optional :: sec
    type(TDate),intent(inout) :: date
    type(TDate) :: new_date
    integer :: cyear,cmonth,cday,chour,cmin
    real(DEFAULT_REAL_BYTE) :: csec
#ifdef DEBUG
    integer :: mday(12)
#endif

    if (present(year)) then
#ifdef DEBUG
    if(year<base_year) then
       write(0,*) 'error: year < base_year'
       stop
    endif
#endif
       cyear = year
    else
       cyear = date%year
    endif
    if (present(month)) then
       cmonth = month
#ifdef DEBUG
    if(month<1 .or. month>12) then
       write(0,*) 'error: month is out of range'
       stop
    endif
#endif
    else
       cmonth = date%month
    endif
    if (present(day)) then
#ifdef DEBUG
       mday=days_month(year)
       if(day<0 .or. day>mday(date%month)) then
          write(0,*) 'error: day is out of range'
          stop
       endif
#endif
       cday = day
    else
       cday = date%day
    endif
    if (present(hour)) then
#ifdef DEBUG
       if(hour<0 .or. hour>23) then
          write(0,*) 'error: hour is out of range'
          stop
       endif
#endif
       chour = hour
    else
       chour = date%hour
    endif
    if (present(min)) then
#ifdef DEBUG
       if(min<0 .or. min>59) then
          write(0,*) 'error: min is out of range'
          stop
       endif
#endif
       cmin = min
    else
       cmin = date%min
    endif
    if (present(sec)) then
#ifdef DEBUG
       if(sec<0 .or. sec>59) then
          write(0,*) 'error: sec is out of range'
          stop
       endif
#endif
       csec = sec
    else
       csec = date%sec
    endif

    new_date = init_date1(cyear,cmonth,cday,chour,cmin,csec)

  end function replace_date


!----------- for TClimdate

! new_date = date1 + date2

  function climdate_plus_climdate(date1,date2) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date1,date2
    type(TClimdate) :: new_date

    new_date%tot_sec = date1%tot_sec + date2%tot_sec

    call update_climdate(new_date)

  end function climdate_plus_climdate

!  new_date = date + timediff
  function climdate_plus_timediff(date,timediff) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    type(TTimediff),intent(in) :: timediff
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec + timediff%diff_sec 

    call update_climdate(new_date)

  end function climdate_plus_timediff

! new_date = timediff + date
  function timediff_plus_climdate(timediff,date) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    type(TTimediff),intent(in) :: timediff
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec + timediff%diff_sec 

    call update_climdate(new_date)

  end function timediff_plus_climdate

! new_date = date - timediff
  function climdate_minus_timediff(date,timediff) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    type(TTimediff),intent(in) :: timediff
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec - timediff%diff_sec 

#ifdef DEBUG
    if (new_date%tot_sec < 0) then
       write(0,*) 'error: tot_sec < 0'
       stop
    endif
#endif

    call update_climdate(new_date)

  end function climdate_minus_timediff


!  timediff = date1 - date2
  function climdate_minus_climdate(date1,date2) result(timediff)
    implicit none
    type(TClimdate),intent(in) :: date1,date2
    type(TTimediff) :: timediff

    timediff%diff_sec = date1%tot_sec - date2%tot_sec

    call update_timediff(timediff)

  end function climdate_minus_climdate
  

! new_date = date * val

  function climdate_times_double(date, val) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    real(8),intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec*val

    call update_climdate(new_date)

  end function climdate_times_double

  function climdate_times_single(date, val) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    real(4),intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec*dble(val)

    call update_climdate(new_date)

  end function climdate_times_single

  function climdate_times_int(date, val) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    integer,intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec*dble(val)

    call update_climdate(new_date)

  end function climdate_times_int


  function double_times_climdate(val, date) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    real(8),intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec*val

    call update_climdate(new_date)

  end function double_times_climdate

  function single_times_climdate(val, date) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    real(4),intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec*dble(val)

    call update_climdate(new_date)

  end function single_times_climdate

  function int_times_climdate(val, date) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    integer,intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec*dble(val)

    call update_climdate(new_date)

  end function int_times_climdate


! new_date = date / val

  function climdate_div_double(date, val) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    real(8),intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec/val

    call update_climdate(new_date)

  end function climdate_div_double


  function climdate_div_single(date, val) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    real(4),intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec/dble(val)

    call update_climdate(new_date)

  end function climdate_div_single


  function climdate_div_int(date, val) result(new_date)
    implicit none
    type(TClimdate),intent(in) :: date
    integer,intent(in) :: val
    type(TClimdate) :: new_date

    new_date%tot_sec = date%tot_sec/dble(val)

    call update_climdate(new_date)

  end function climdate_div_int

  
! replace date member

  function replace_climdate(date,month,day,hour,min,sec,sum_year) &
       result(new_date)
    implicit none
    integer,intent(in),optional :: sum_year,month,day,hour,min
    real(DEFAULT_REAL_BYTE),intent(in),optional :: sec
    type(TClimdate),intent(in) :: date
    type(TClimdate) :: new_date
    integer :: cyear,cmonth,cday,chour,cmin
    real(DEFAULT_REAL_BYTE) :: csec
#ifdef DEBUG
    integer :: mday(12)
#endif

    if (present(sum_year)) then
       cyear = sum_year
    else
       cyear = date%sum_year
    endif
    if (present(month)) then
#ifdef DEBUG
    if(month<1 .or. month>12) then
       write(0,*) 'error: month is out of range'
       stop
    endif
#endif
       cmonth = month
    else
       cmonth = date%month
    endif
    if (present(day)) then
#ifdef DEBUG
       mday=days_month_clim()
       if(day<0 .or. day>mday(date%month)) then
          write(0,*) 'error: day is out of range'
          stop
       endif
#endif
       cday = day
    else
       cday = date%day
    endif
    if (present(hour)) then
#ifdef DEBUG
       if(hour<0 .or. hour>23) then
          write(0,*) 'error: hour is out of range'
          stop
       endif
#endif
       chour = hour
    else
       chour = date%hour
    endif
    if (present(min)) then
#ifdef DEBUG
       if(min<0 .or. min>59) then
          write(0,*) 'error: min is out of range'
          stop
       endif
#endif
       cmin = min
    else
       cmin = date%min
    endif
    if (present(sec)) then
#ifdef DEBUG
       if(sec<0 .or. sec>59) then
          write(0,*) 'error: sec is out of range'
          stop
       endif
#endif
       csec = sec
    else
       csec = date%sec
    endif

    new_date = init_climdate1(cmonth,cday,chour,cmin,csec,cyear)

  end function replace_climdate


!-----------------------  subordinate functions  ---------------------
  subroutine wdhms_from_sec(insec,week,day,hour,min,sec)
    implicit none
    integer,intent(out) :: hour,min,day,week
    real(DEFAULT_REAL_BYTE),intent(out) :: sec
    real(8),intent(in) :: insec
    real(8) :: cday

    cday = dble(int(insec/dble(sec_a_day)))

    week = cday/days_a_week
    day = mod(int(cday),days_a_week)
    sec = insec-cday*dble(sec_a_day)
    min = mod(int(sec/sec_a_min), min_an_hour)
    hour = sec/sec_an_hour

  end subroutine wdhms_from_sec

  function cal_tot_sec(date) result(tot_sec)
    implicit none
    type(TDate),intent(in) :: date
    integer :: cyear
    real(8) :: tot_sec,sday

    cyear = base_year
    sday = 0.d0
    do cyear = base_year, date%year - 1
       sday = sday + dble(days_year(cyear))
    end do

    sday = sday + dble(year_day(date)) - 1.d0

    tot_sec = sday*dble(sec_a_day) + dble(date%hour*sec_an_hour) &
         + dble(date%min*sec_a_min) + dble(date%sec)

  end function cal_tot_sec

  function cal_tot_sec_clim(date) result(tot_sec)
    implicit none
    type(TClimdate) :: date
    real(8) :: tot_sec,sday
    integer :: cmonth

    sday = dble(date%sum_year)*dble(days_a_year_clim)
    sday = sday + dble(year_day_clim(date)) - 1.d0


    tot_sec = sday*dble(sec_a_day) + dble(date%hour*sec_an_hour) &
         + dble(date%min*sec_a_min) + dble(date%sec)

  end function cal_tot_sec_clim

  function year_day(date)
    implicit none
    type(TDate),intent(in) :: date
    integer :: year_day,mday(12)

    if(date%month > 1) then
       mday = days_month(date%year)
       year_day = sum(mday(1:date%month-1))
       year_day = year_day + date%day
    else
       year_day = date%day
    end if

  end function year_day

  function year_day_clim(date) result(year_day)
    implicit none
    type(TClimdate),intent(in) :: date
    integer :: year_day,mday(12)

    if(date%month > 1) then
       mday = days_month_clim()
       year_day = sum(mday(1:date%month-1))
       year_day = year_day + date%day
    else
       year_day = date%day
    end if


  end function year_day_clim

  function year_sec(date)
    implicit none
    type(TDate),intent(in) :: date
    real(8) :: year_sec
    
    year_sec = dble(year_day(date)-1)*dble(sec_a_day)&
         + dble(date%hour)*dble(sec_an_hour) +&
         dble(date%min)*dble(sec_a_min) + dble(date%sec)

  end function year_sec

  function year_sec_clim(date) result(year_sec)
    implicit none
    type(TClimdate),intent(in) :: date
    real(8) :: year_sec
    
    year_sec = dble(year_day_clim(date)-1)*dble(sec_a_day)&
         + dble(date%hour)*dble(sec_an_hour) +&
         dble(date%min)*dble(sec_a_min) + dble(date%sec)

  end function year_sec_clim

  function leap_year(year)
    implicit none
    integer,intent(in) :: year
    logical :: leap_year

    if ((mod(year,4) == 0 .and. mod(year,100) /= 0) &
         .or. mod(year,1000) == 0) then
       leap_year = .true.
    else
       leap_year = .false.
    end if
    
  end function leap_year

  function  year_day_to_date(year,year_day) result(date)
    implicit none
    integer,intent(in) :: year,year_day
    type(TDate) :: date
    integer :: mday(12),cyear,cmonth,cday,sum_day

    cyear = year
    cday = year_day
    do
       if (cday < 1) then
          cyear = cyear - 1
#ifdef DEBUG
          if(cday < base_year) then
             write(0,*) 'error: year < base_year'
             stop
          endif
#endif
          cday = cday + days_year(cyear)
       elseif (cday > days_year(cyear)) then
          cday = cday - days_year(cyear)
          cyear = cyear + 1
       else
          exit
       end if
    end do

    mday=days_month(cyear)
    date%year = cyear
   
    if (cday <= mday(1)) then
       date%day = cday
       date%month = 1
       return
    end if
    
    sum_day = mday(1)
    cmonth = 1
    do while(cday > sum_day)
       cmonth = cmonth + 1
       sum_day = sum_day + mday(cmonth)
    end do

    date%month = cmonth
    date%day = cday-sum(mday(1:cmonth-1))

!    date%tot_sec = cal_tot_sec(date)

  end function year_day_to_date

  function year_day_to_climdate(year_day,sum_year) result(date)
    implicit none
    integer,intent(in) :: year_day,sum_year
    type(TClimdate) :: date
    integer :: cday,cyear,cmonth,mday(12),sum_day

    cyear = sum_year
    cday = year_day
    do
       if (cday < 1) then
          cyear = cyear - 1
          cday = cday + days_a_year_clim
       elseif (cday > days_a_year_clim) then
          cyear = cyear + 1
          cday = cday - days_a_year_clim
       else
          exit
       end if
    end do

    mday=days_month_clim()
    date%sum_year = cyear
   
    if (cday <= mday(1)) then
       date%day = cday
       date%month = 1
       return
    end if
    
    sum_day = mday(1)
    cmonth = 1
    do while(cday > sum_day)
       cmonth = cmonth + 1
       sum_day = sum_day + mday(cmonth)
    end do

    date%month = cmonth
    date%day = cday-sum(mday(1:cmonth-1))

  end function year_day_to_climdate

  function days_month(year)
    implicit none
    integer,intent(in) :: year
    integer :: days_month(12)

    days_month(1:12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    if (leap_year(year)) then
       days_month(2) = 29
    else
       days_month(2) = 28
    end if

  end function days_month
  
  function days_month_clim()
    implicit none 
    integer :: days_month_clim(12)

    select case(days_a_year_clim)
    case(360)
       days_month_clim(1:12) = 30
    case(365)
       days_month_clim(1:12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    case(366)
       days_month_clim(1:12) = (/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
#ifdef DEBUG
    case default
       write(0,*) 'error: days_a_year_clim must be 360, 365 or 366'
       stop
#endif
    end select

  end function days_month_clim
    
  function days_year(year)
    integer,intent(in) :: year
    integer :: days_year

    if (leap_year(year)) then
       days_year = 366
    else
       days_year = 365
    end if

  end function days_year

  function year_sec_to_date(year,year_sec) result(date)
    implicit none
    integer,intent(in) :: year
    real(8),intent(in) :: year_sec
    real(DEFAULT_REAL_BYTE) :: sec_the_day
    type(TDate) :: date
    integer :: iyear_sec
    iyear_sec = int(year_sec)

    sec_the_day = mod(iyear_sec, sec_a_day)+year_sec-iyear_sec

    date = year_day_to_date(year,int(year_sec)/sec_a_day+1)

    date%hour = sec_the_day / sec_an_hour 
    date%min = sec_the_day/sec_a_min - date%hour*min_an_hour
    date%sec = sec_the_day - dble(date%min*sec_a_min) -&
         dble(date%hour*sec_an_hour)
    date%tot_sec = cal_tot_sec(date)

  end function year_sec_to_date

  function year_sec_to_climdate(year_sec,sum_year) result(date)
    implicit none
    integer,intent(in) :: sum_year
    real(8),intent(in) :: year_sec
    real(DEFAULT_REAL_BYTE) :: sec_the_day
    type(TClimdate) :: date
    integer :: iyear_sec
    iyear_sec = int(year_sec)

    sec_the_day = mod(iyear_sec, sec_a_day)+year_sec-iyear_sec

    date = year_day_to_climdate(int(year_sec)/sec_a_day+1,sum_year)

    date%hour = sec_the_day / sec_an_hour 
    date%min = sec_the_day/sec_a_min - date%hour*min_an_hour
    date%sec = sec_the_day - dble(date%min*sec_a_min) -&
         dble(date%hour*sec_an_hour)
    date%tot_sec = cal_tot_sec_clim(date)

  end function year_sec_to_climdate

#ifndef NO_MPI
     subroutine bcast_tdate(date_val)
     
     use mod_mpi
     implicit none
     type(Tdate) date_val
     
     integer idum(6)
     real(DEFAULT_REAL_BYTE) :: sc

     if(ip==imaster) then
       idum(1) = date_val%year
       idum(2) = date_val%month
       ! idum(3) = date_val%hour !BUG
       ! idum(4) = date_val%day  !BUG
       idum(3) = date_val%day   !fix by LQH
       idum(4) = date_val%hour  !fix by LQH
       idum(5) = date_val%min
       sc = date_val%sec
      endif
      
      call bcast_int(idum,5)
#if (DEFAULT_REAL_BYTE == 8) 
      call bcast_dble1(sc,1)
#elif (DEFAULT_REAL_BYTE == 4)
      call bcast_real1(sc,1)
#else
      write(*,*) 'DEFAULT_REAL_BYTE is bad'
#endif

     if(ip/=imaster) then
       date_val=init_date1(idum(1),idum(2),idum(3),idum(4),idum(5),sc)
     endif
     
     end subroutine bcast_tdate
#endif

end module datetime

#ifdef TEST
program test
 use datetime
 implicit none
 type(TDate) :: date1,date2,date
 type(TTimediff) :: dt,dt2
 type(TClimdate) :: cdate1,cdate2
 real :: t,t1,t2,coef
 character(len=3) :: month_label(12) = (/'Jan','Feb','Mar','Apr','May',&
       'Jun','Jul','Aug','Sep','Oct','Nov','Dec' /)

! initialize TDate
 date1 = init_date1(2007,9,15,10)
 date2 = init_date1(2008,3,11,22,31,11.)


! initialize TTimediff
 dt = timediff(day=12,hour=10,min=1,sec=4.1)
 dt2= date2 - date1

! (TTimediff = TDate - TDate)
! (Do NOT dt2%diff_sec = date2%tot_sec - date1%tot_sec )

 t1 = 23.
 t2 = 11.

 print *,'Stating date'
 call print_date(date1)
 print *,'Value =',t1
 print *,''
 print *,'Last date'
 call print_date(date2)
 print *,'Value =',t2
 print *,''
 print *,'Time step = ',dt%diff_sec,'sec'
 print *,''
 print *,'Total time =',dt2%diff_sec,'sec'
 date = date1
! when comparing dates, refer member "tot_sec"
 do while (date%tot_sec < date2%tot_sec)

    coef = (date2 - date)/dt2
! (REAL = TTimediff / TTimediff)
! (identical : coef = (data2%tot_sec - data%tot_sec)/dt2%diff_sec )

    t = t2 - coef*(t2 - t1)

    call print_date(date)
    print *,'Value =',t

    date = date + dt   
! (TDate = TDate + TTimediff)
 end do

 contains 
   subroutine print_date(date)
     type(TDate) :: date

     write(*,fmt="(i4.4,1x,a,1x,i2,1x,i2.2,':',i2.2,':',f4.1)")&
          date%year,month_label(date%month),date%day,date%hour,&
          date%min,date%sec

   end subroutine print_date
end program test
#endif
