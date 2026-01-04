c**** NCF_SET_TSTEP
c
      subroutine ncf_set_tstep(begin_date,begin_time,
     &                                       ending_date,ending_time)
      use ncf_iomod
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the global file attributes for the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         begin_date  I model begin date (YYJJJ)
c         begin_time  R model begin time
c         ending_date I model end date (YYJJJ)
c         ending_time R model end time
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.inc'
      include 'namelist.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer      begin_date
      real         begin_time
      integer      ending_date
      real         ending_time
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer this_date, date_now, num_steps, i
      real    time_now, this_hour, this_min
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c  --- initialize to start of simulation ---
c
c
c  --- initialize to start of simulation ---
c
      ncf_nsteps = 0
      date_now = begin_date
      time_now = begin_time
c
c  --- walk through time as the the simulation would output,
c      this is to figure out how many steps for allocation ---
c
 111  continue
      this_hour = AINT(time_now/100.)
      this_min = AMOD(time_now,100.)
      time_now = 100.*this_hour + this_min
      this_date = date_now
      if(time_now .GE. 2400.) then
        time_now = time_now - 2400.
        this_date = this_date + 1
        if( MOD(this_date,1000) .GT. 365 ) then
           if( MOD(INT(this_date/1000),4) .EQ. 0 ) then
              if( MOD(this_date,1000) .EQ. 367 )
     &                   this_date = (INT(this_date/1000)+1)*1000 + 1
           else
              this_date = (INT(this_date/1000)+1)*1000 + 1
           endif
        endif
      endif
      ncf_nsteps = ncf_nsteps + 1
      call uptime(time_now,date_now,dtout*60.)
      if(this_date .LT. ending_date) goto 111
      if(this_date .EQ. ending_date .AND. 
     &                     time_now .LT. ending_time - 0.01) goto 111
c
c  --- adjust for end of day ---
c
      if( time_now .GT. ending_time ) ncf_nsteps = ncf_nsteps - 1
c
c  --- call routine to allocate allocate the array ----
c
      call ncf_alloc_tstep(ncf_nsteps)
c
c  --- now walk through time again to fill the NetCDF variable ---
c
      num_steps = 0
      date_now = begin_date
      time_now = begin_time
c
c  --- walk through time as the the simulation would output,
c      this is to figure out how many steps for allocation ---
c
 222  continue
      this_hour = AINT(time_now/100.)
      this_min = AMOD(time_now,100.)
      time_now = 100.*this_hour + this_min
      this_date = date_now
      if(time_now .GE. 2400.) then
        time_now = time_now - 2400.
        this_date = this_date + 1
        if( MOD(this_date,1000) .GT. 365 ) then
           if( MOD(INT(this_date/1000),4) .EQ. 0 ) then
              if( MOD(this_date,1000) .EQ. 367 )
     &                   this_date = (INT(this_date/1000)+1)*1000 + 1
           else
              this_date = (INT(this_date/1000)+1)*1000 + 1
           endif
        endif
      endif
      num_steps = num_steps + 1
      if( num_steps .GT. ncf_nsteps ) goto 333
      ncf_tflag(1,num_steps) = 100000*(Start_Date_Hour(1)/100) +
     &                         this_date     
      ncf_tflag(2,num_steps) = time_now*100
      if( num_steps .GT. 1 ) then
         ncf_etflag(1,num_steps-1) = 100000*(Start_Date_Hour(1)/100) +
     &                               this_date
         ncf_etflag(2,num_steps-1) = time_now*100
      endif
      call uptime(time_now,date_now,dtout*60.)
      goto 222
c
c  --- reduce count by 1 because the ending is not a start date/time ---
c
 333  continue
c
c  --- finish with the ending time
c
      if( ncf_nsteps .GT. 0 ) then
         ncf_etflag(1,ncf_nsteps) = 100000*(Start_Date_Hour(1)/100) + ending_date
         ncf_etflag(2,ncf_nsteps) = ending_time*100
      endif
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
