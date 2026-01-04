C**** NCF_READINP_3D
c
      subroutine ncf_readinp_3d(igrd,ncol,nrow,nlay,height,phptim,hnext,
     &                             press,ppptim,pnext,windu,puptim,unext,
     &                                    windv,pvptim,vnext,tempk,ptptim,
     &                                            tnext,water,pwptim,wnext)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routines reads the next hour of data in the 3D met file and
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
c        height              layer interface height field (m)
c        press               layer midpoint pressure field (mb)
c        windu               u-component wind field (m/s)
c        windv               u-component wind field (m/s)
c        tempk               temperature field (K)
c        water               water vapor field (ppm)
c
c     Output arguments:
c        phptim              time-rate change of layer interface height (m/s)
c        hnext                next layer interface height field (m)
c        ppptim              time-rate change of pressure (mb/s)
c        pnext                next pressure field (mb)
c        puptim              time-rate change of u-component wind (m/s2)
c        unext                next u-component wind field (m/s)
c        pvptim              time-rate change of v-component wind (m/s2)
c        vnext                next v-component wind field (m/s)
c        psptim              time-rate change of surface temperature (K/s)
c        tsnext               next surface temperature field (m/s)
c        ptptim              time-rate change of 3-D temperature (K/s)
c        tnext                next 3-D temperature field (m/s)
c        pwptim              time-rate change of 3-D water vapor (ppm/s)
c        wnext                next 3-D water vapor field (ppm)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     07/23/19   Staggered wind flag is now grid-specific
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
      real    height(ncol,nrow,nlay)
      real    phptim(ncol,nrow,nlay)
      real    hnext(ncol,nrow,nlay)
      real    press(ncol,nrow,nlay)
      real    ppptim(ncol,nrow,nlay)
      real    pnext(ncol,nrow,nlay)
      real    windu(ncol,nrow,nlay)
      real    puptim(ncol,nrow,nlay)
      real    unext(ncol,nrow,nlay)
      real    windv(ncol,nrow,nlay)
      real    pvptim(ncol,nrow,nlay)
      real    vnext(ncol,nrow,nlay)
      real    tempk(ncol,nrow,nlay)
      real    ptptim(ncol,nrow,nlay)
      real    tnext(ncol,nrow,nlay)
      real    water(ncol,nrow,nlay)
      real    pwptim(ncol,nrow,nlay)
      real    wnext(ncol,nrow,nlay)
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
      integer       hdate, this_tstep, i, j, k
      real          whr, wmn, atim, htim
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if( .NOT. is_netcdf_i3dmet(igrd) ) goto 9999
c
      write(action,'(A,I3)') 'Reading 3D met file for grid ',igrd
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
      this_tstep = ncf_get_tstep(i3dmet(igrd),action,hdate,htim,
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
      this_var = "z"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,hnext)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- pressure ---
c
      this_var = "pressure"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,pnext)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- temperature ---
c
      this_var = "temperature"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,tnext)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- humidity ---
c
      this_var = "humidity"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,wnext)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- U-component of wind ---
c
      this_var = "uwind"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,unext)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- V-component of wind ---
c
      this_var = "vwind"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,vnext)
      if( ierr .NE. NF_NOERR) goto 7001
c
      if (.not.lstagw(igrd)) call cvtwind(ncol,nrow,nlay,unext,vnext)
c
c-----Calculate time rates of change for 3D met fields
c
      call timrates(ncol,nrow,nlay,height,hnext,phptim)
      call timrates(ncol,nrow,nlay,press,pnext,ppptim)
      call timrates(ncol,nrow,nlay,tempk,tnext,ptptim)
      call timrates(ncol,nrow,nlay,water,wnext,pwptim)
      call timrates(ncol,nrow,nlay,windu,unext,puptim)
      call timrates(ncol,nrow,nlay,windv,vnext,pvptim)
c
c  --- write message to out file ---
c
      write(iout,'(a40,f7.0,i8.5,a,i3)') 'Read 3D met file at ',htim,hdate,
     &                                                         ' grid',igrd
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
      write(iout,'(//,a)') 'ERROR in NCF_READINP_3D:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_READINP_3D:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      return
      end
