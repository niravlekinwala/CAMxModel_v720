C**** NCF_READINP_KV
c
      subroutine ncf_readinp_kv(igrd,ncol,nrow,nlay,rkv,pkptim,knext)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routines reads the next hour of data in the KV met file and
c   calculates the time/rate change based on the data for the current 
c   hour.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        rkv                 current vertical exchange coefficient (m2/s)
c
c     Output arguments:
c        pkptim              time-rate change of 3-D vertical diff (m2/s2)
c        knext               next 3-D vertical diffusivity field (m2/s)
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
      include 'camx.inc'
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrd
      integer ncol
      integer nrow
      integer nlay
      real    rkv(ncol,nrow,nlay)
      real    pkptim(ncol,nrow,nlay)
      real    knext(ncol,nrow,nlay)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
      integer ncf_get_tstep
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      character*20  this_var
      integer       this_varid, ierr, data_start(4), data_count(4)
      integer       hdate, this_tstep, i, j, k
      integer       this_time_tflag, this_time_etflag
      real          whr, wmn, atim, htim
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if( .NOT. is_netcdf_ikv(igrd) ) goto 9999
c
      write(action,'(A,I3)') 'Reading KV met file for grid ',igrd
c
c --- Find the date/time at the next update interval ---
c
      whr = aint(time/100.)
      wmn = anint(amod(time,100.))
      atim = 100.*whr + wmn
      htim = 100.*(whr + aint((wmn + dtinp)/60.)) +
     &             amod((wmn + dtinp),60.)
      hdate = date
      if (htim.ge.2400.) then
        htim = anint(htim - 2400.)
        hdate = hdate + 1
        if( MOD(hdate,1000) .GT. 365 ) then
           if( MOD(INT(hdate/1000),4) .EQ. 0 ) then
              if( MOD(hdate,1000) .EQ. 367 )
     &           hdate = (INT(hdate/1000)+1)*1000 + 1
           else
              hdate = (INT(hdate/1000)+1)*1000 + 1
           endif
        endif
      endif
c
c   --- get the timestep that encompasses this date/time ---
c
      this_tstep = ncf_get_tstep(ikv(igrd),action,hdate,htim,
     &               this_time_tflag,this_time_etflag,.FALSE.,.FALSE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol
      data_start(2) = 1
      data_count(2) = nrow
      data_start(3) = 1
      data_count(3) = nlay
      data_start(4) = this_tstep
      data_count(4) = 1
c
c ---- layer heights ---
c
      this_var = "kv"
      ierr = nf_inq_varid(ikv(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(ikv(igrd),this_varid,data_start,
     &                                          data_count,knext)
      if( ierr .NE. NF_NOERR) goto 7001
c
c-----Calculate time rates of change for KV met fields
c
      call timrates(ncol,nrow,nlay,rkv,knext,pkptim)
c
c  --- write message to out file ---
c
      write(iout,'(a40,f7.0,i8.5,a,i3)') 
     &    'Read KV met file at ',htim,hdate,' grid',igrd
      call flush(iout)
c
c  --- successful completion ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_READINP_KV:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_READINP_KV:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      return
      end
