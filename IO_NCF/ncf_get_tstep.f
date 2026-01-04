c**** NCF_GET_TSTEP
c
      function ncf_get_tstep(iounit,action,this_date,this_time,
     &         this_time_tflag,this_time_etflag,ignore_date,lstrict)
      use filunit
      implicit none
      integer ncf_get_tstep
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine retrieves the time step flags from the NetCDF file
c   and looks for a time span that encompasses the current data/time.
c   It returns the index in the time flag array.
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         iounit            I NCF file ID
c         action            C string that describes file being read
c         this_date         I current model date (YYJJJ)
c         this_time         R current model time
c         this_time_tflag   R current model time
c         this_time_etflag  R current model time
c         ignore_date       L .FALSE. if not checking date
c         lstrict           L .TRUE. if must be strictly less than endtime
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
      integer       this_date
      real          this_time
      real          this_time_tflag
      real          this_time_etflag
      logical       ignore_date
      logical       lstrict
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
      character*10 this_var
      integer      tflag_varid, etflag_varid, this_tstep, tstep_in
      integer      data_start(3), data_count(3), nvars_in, ierr, i
      integer      tflag_first(2), tflag_last(2)
      integer      etflag_first(2), etflag_last(2)
      real         date_time, date_time_tflag, date_time_etflag
      logical      lread_ok
c
      integer,     allocatable, dimension(:,:) :: tflag_in
      integer,     allocatable, dimension(:,:) :: etflag_in
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
      ncf_get_tstep = -9
c
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
      date_time = REAL(this_date)+this_time/2400.
c
c  --- get each time flag pair and check that 
c
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
         date_time_tflag = REAL(this_date)+REAL(tflag_in(2,1))/240000.
         date_time_etflag = REAL(this_date)+REAL(etflag_in(2,1))/240000.
         if( etflag_in(2,1) .EQ. 0. ) date_time_etflag = date_time_etflag + 1.
      else
         if( MOD(tflag_in(1,1)/1000,100) .LE. 80 ) then
            date_time_tflag = REAL(tflag_in(1,1)-2000000)+REAL(tflag_in(2,1))/240000.
            date_time_etflag = REAL(etflag_in(1,1)-2000000)+REAL(etflag_in(2,1))/240000.
         else
            date_time_tflag = REAL(tflag_in(1,1)-1900000)+REAL(tflag_in(2,1))/240000.
            date_time_etflag = REAL(etflag_in(1,1)-1900000)+REAL(etflag_in(2,1))/240000.
         endif
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
      if( lstrict ) then
          if( date_time .GE. date_time_tflag .AND. 
     &                        date_time .LT. date_time_etflag ) then
             ncf_get_tstep = this_tstep
             this_time_tflag = date_time_tflag
             this_time_etflag = date_time_etflag
             goto 9999
          endif
      else
          if( date_time .GE. date_time_tflag .AND. 
     &                        date_time .LE. date_time_etflag ) then
             ncf_get_tstep = this_tstep
             this_time_tflag = date_time_tflag
             this_time_etflag = date_time_etflag
             goto 9999
          endif
      endif
      goto 111
c
 222  continue
      write(iout,'(//,a)') 'ERROR in NCF_GET_TSTEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Did not find suitable timesteps in the file ',
     &                                     'for this time period.'
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
      write(iout,'(//,a)') 'ERROR in NCF_GET_TSTEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_GET_TSTEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get global attributes for time span.',
     &                                                ' Looking for: '
      write(iout,'(10X,A)') 'TSTEP'
      write(iout,'(10X,A)') 'NVARS'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
