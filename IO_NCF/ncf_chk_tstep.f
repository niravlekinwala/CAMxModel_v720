c**** NCF_CHK_TSTEP
c
      subroutine ncf_chk_tstep(iounit,action,begin_date,begin_time,
     &                           ending_date,ending_time,ignore_date)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine checks that the time span in the file encompasses
c   the simulation period and that there is gridded data for each
c   tome period that will be accessed.
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         iounit      I NCF file ID
c         action      C string that describes file being read
c         begin_date  I model begin date (YYJJJ)
c         begin_time  R model begin time
c         ending_date I model end date (YYJJJ)
c         ending_time R model end time
c         ignore_date L .FALSE. if date should not be checked
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
      include 'camx.prm'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer       iounit
      character*(*) action
      integer       begin_date
      real          begin_time
      integer       ending_date
      real          ending_time
      logical       ignore_date
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 this_var
      integer       sdate_mod, sdate_in, stime_mod, stime_in, nsteps
      integer       tstep_in, nvars_in, tflag_varid, etflag_varid
      integer       data_start(3), data_count(3), this_tstep, ierr
      integer       tflag_first(2), tflag_last(2)
      integer       etflag_first(2), etflag_last(2)
      real          date_time_beg, date_time_end
      real          date_time_tflag, date_time_etflag
      logical       lread_ok, lcheck_begtime, lcheck_endtime
c
      integer,      allocatable, dimension(:,:) :: tflag_in
      integer,      allocatable, dimension(:,:) :: etflag_in
c
c-----------------------------------------------------------------------
c    Parameters:
c-----------------------------------------------------------------------
c
      real FUZZ
      parameter( FUZZ = 0.01 )
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- assume all variables are available ----
c
      lread_ok = .TRUE.
c
c  --- get the time variables needed to check ---
c
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'SDATE', sdate_in)
      if( ierr .NE. NF_NOERR ) lread_ok = .FALSE.
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'STIME', stime_in)
      if( ierr .NE. NF_NOERR ) lread_ok = .FALSE.
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'TSTEP', tstep_in)
      if( ierr .NE. NF_NOERR ) lread_ok = .FALSE.
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NVARS', nvars_in)
      if( ierr .NE. NF_NOERR ) lread_ok = .FALSE.
      if( .NOT. lread_ok ) goto 7001
c
c  --- allocate the arrays for the time flags ---
c
      allocate( tflag_in(2,nvars_in))
      allocate( etflag_in(2,nvars_in))
c
c  --- calculate date/time together ---
c
      if( ignore_date ) then
         if( begin_date/1000 .LE. 80 ) then
            date_time_beg = REAL(2000000+begin_date)+begin_time/2400.
            date_time_end = REAL(2000000+begin_date)+ending_time/2400.
         else
            date_time_beg = REAL(1900000+begin_date)+begin_time/2400.
            date_time_end = REAL(1900000+ending_date)+ending_time/2400.
         endif
      else
         if( begin_date/1000 .LE. 80 ) then
             date_time_beg = REAL(2000000+begin_date)+begin_time/2400.
             date_time_end = REAL(2000000+ending_date)+ending_time/2400.
          else
             date_time_beg = REAL(1900000+begin_date)+begin_time/2400.
             date_time_end = REAL(1900000+ending_date)+ending_time/2400.
          endif
      endif
c
c  --- get each time flag pair and check that 
c
      lcheck_begtime = .FALSE.
      lcheck_endtime = .FALSE.
      this_var = 'TFLAG'
      ierr = nf_inq_varid(iounit,this_var,tflag_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      this_var = 'ETFLAG'
      ierr = nf_inq_varid(iounit,this_var,etflag_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      data_start(1) = 1
      data_count(1) = 2
      data_start(2) = 1
      data_count(2) = nvars_in
      this_tstep = 0
 111  continue
      this_tstep = this_tstep + 1
      data_start(3) = this_tstep
      data_count(3) = 1
      ierr = nf_get_vara_int(iounit,tflag_varid,data_start,data_count,tflag_in)
      if( ierr .NE. NF_NOERR ) goto 222
      ierr = nf_get_vara_int(iounit,etflag_varid,data_start,data_count,etflag_in)
      if( ierr .NE. NF_NOERR ) goto 222
c
c  --- calculate date/time together ---
c
      if( ignore_date ) then
         if( begin_date/1000 .LE. 80 ) then
            date_time_tflag = REAL(2000000+begin_date)+REAL(tflag_in(2,1))/240000.
            date_time_etflag = REAL(2000000+begin_date)+REAL(etflag_in(2,1))/240000.
         else
            date_time_tflag = REAL(1900000+begin_date)+REAL(tflag_in(2,1))/240000.
            date_time_etflag = REAL(1900000+begin_date)+REAL(etflag_in(2,1))/240000.
         endif
      else
         date_time_tflag = REAL(tflag_in(1,1))+REAL(tflag_in(2,1))/240000.
         date_time_etflag = REAL(etflag_in(1,1))+REAL(etflag_in(2,1))/240000.
      endif
      if( this_tstep .EQ. 1 ) then
          tflag_first(1) = tflag_in(1,1)
          tflag_first(2) = tflag_in(2,1)
          etflag_first(1) = etflag_in(1,1)
          etflag_first(2) = etflag_in(2,1)
      endif
      tflag_last(1) = tflag_in(1,1)
      tflag_last(2) = tflag_in(2,1)
      etflag_last(1) = etflag_in(1,1)
      etflag_last(2) = etflag_in(2,1)
c
c  --- check if this time period covers the beginning episode time ---
c
      if( date_time_beg .GE. date_time_tflag .AND. date_time_beg .LE.
     &                              date_time_etflag) lcheck_begtime = .TRUE.
c
c  --- check if this time period covers the ending episode time ---
c
      if( date_time_end .GE. date_time_tflag .AND. date_time_end .LE.
     &                              date_time_etflag) lcheck_endtime = .TRUE.
c
      if( lcheck_endtime .AND. lcheck_endtime ) goto 9999
      goto 111
c
 222  continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_TSTEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Did not find suitable timesteps in the file ',
     &                                     'for this simulation period.'
      write(iout,'(A,4I10)') 'First timestep: ',tflag_first(1), tflag_first(2),
     &                                        etflag_first(1), etflag_first(2)
      write(iout,'(A,4I10)') 'Last timestep : ',tflag_last(1), tflag_last(2),
     &                                        etflag_last(1), etflag_last(2)
      write(iout,'(2A)') 'Check that the file contains all hours ',
     &                                               'to be simulated.'
      call camxerr()
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_TSTEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_TSTEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get global attributes for time span.',
     &                                                ' Looking for: '
      write(iout,'(10X,A)') 'SDATE'
      write(iout,'(10X,A)') 'STIME'
      write(iout,'(10X,A)') 'TSTEP'
      write(iout,'(10X,A)') 'NVARS'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
