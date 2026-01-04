C**** GETIME.F 
c
      subroutine getime( this_date, this_time )
c
c-----------------------------------------------------------------------
c
c   Description:
c
c     This routine returns the current date/time as integer values in 
c     IOAPI format.
c
c   Arguments:
c
c     Outputs:
c       this_date    I   date (YYYYJJJ)
c       this_time    I   time (HHMMSS)
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   Argument declaration:
c-----------------------------------------------------------------------
c
      integer this_date
      integer this_time
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*8  cdate
      character*10 ctime
      character*5  czone
      integer*4    imon, iday, iyear, ihr, imin, isec
      integer*4    iarray(8)
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c   ---- call routine to get date/time as integers ---
c
      call date_and_time(cdate, ctime, czone, iarray) 
      iyear = iarray(1)
      imon = iarray(2)
      iday = iarray(3)
      ihr = iarray(5)
      imin = iarray(6)
      isec = iarray(7)
      iyear = MOD(iyear,100)
      this_date = iyear*10000 + imon*100 + iday
      call juldate(this_date)
      this_date = 2000000 + iyear*1000 + MOD(this_date,1000)
      this_time = ihr*10000 + imin*100 + isec
c
      goto 9999
c
c-----------------------------------------------------------------------
c   Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
