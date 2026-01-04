**** NCF_METINIT_2D
c
      subroutine ncf_metinit_2d(igrd,ncol,nrow,tsurf,snow,snowage,snowalb)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routines initalizes the met fields by reading the NetCDF
c   met files and loading the data into the global arrays.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c
c     Output arguments:
c        tsurf               surface temperature field (K)
c        snow                snow cover water equivalent (m)
c        snowage             snow cover age (hr)
c        snowalb             snow cover albedo (unitless)
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
      real    tsurf(ncol,nrow)
      real    snow(ncol,nrow)
      real    snowage(ncol,nrow)
      real    snowalb(ncol,nrow)
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
      integer       this_time_tflag, this_time_etflag
      integer       this_tstep, i, j, k
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if( .NOT. is_netcdf_i2dmet(igrd) ) goto 9999
c
c  --- initialize fields ---
c
      snow = 0.
      snowage = 0.
      snowalb = 0.
c
      write(action,'(A,I3)') 'Reading 2D met file for grid ',igrd
c
c   --- get the timestep that encompasses this date/time ---
c
      this_tstep = ncf_get_tstep(i2dmet(igrd),action,begdate,begtim,
     &                  this_time_tflag,this_time_etflag,.FALSE.,.FALSE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol
      data_start(2) = 1
      data_count(2) = nrow
      data_start(3) = 1
      data_count(3) = 1
      data_start(4) = this_tstep
      data_count(4) = 1
c
c ---- surface temperature ---
c
      this_var = "sfctemperature"
      ierr = nf_inq_varid(i2dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i2dmet(igrd),this_varid,data_start,
     &                                          data_count,tsurf)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- snow cover depth ---
c
      this_var = "snowewd"
      ierr = nf_inq_varid(i2dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         ierr = nf_get_vara_real(i2dmet(igrd),this_varid,data_start,
     &                                                data_count,snow)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- snow cover age ---
c
      this_var = "snowage"
      ierr = nf_inq_varid(i2dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         ierr = nf_get_vara_real(i2dmet(igrd),this_varid,data_start,
     &                                             data_count,snowage)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
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
      write(iout,'(//,a)') 'ERROR in NCF_METINIT_2D:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_METINIT_2D:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      return
      end
