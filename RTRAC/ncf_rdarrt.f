c***** NCF_RDARRT.F
c
      subroutine ncf_rdarrt(igrid,ndate,ttime,nox,noy,nlay_ems,nrtarsp,emisrt)
      use tracer
      use filunit
      use grid
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions for the RTRAC process
c   and fills the approproate arrays.  The emissions file for one grid
c   but each emissions groups is read. This version is for NetCDF files.
c    Argument descriptions:
c     Outputs:
c       emisrt    R    array of emissions for RTRAC speices
c     Inputs:
c       igrid     I    grid number
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       nox       I    number of columns in grid
c       noy       I    number of rows in grid
c       nlay_ems  I    number of layers in emissions files
c       nrtarsp   I    number of RTRAC species
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c       06/09/18 -gwilson-  original development
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.inc'
      include 'flags.inc'
      include 'netcdf.inc'
c       
c-----------------------------------------------------------------------
c   External functions:
c-----------------------------------------------------------------------
c
      integer istrln
      integer ncf_get_tstep
c       
c-----------------------------------------------------------------------
c   Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      integer ndate
      real    ttime
      integer nox
      integer noy
      integer nlay_ems
      integer nrtarsp
      real    emisrt(nox,noy,nlay_ems,nrtarsp)
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname, action
      character*10  this_var
      integer      nedge, iounit, ierr,  itmp, jtmp
      integer      nlays_in, ispc, data_start(4), data_count(4), this_date
      integer      this_time_tflag, this_time_etflag, this_tstep
      integer      iddm, i, j, k, imap,  this_varid, buffer_offset
      real         emstmp, this_time
      logical      luse
c
      real, allocatable, dimension(:,:,:) :: emsgrd
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
      else
        nedge = 1
      endif
c
c   --- skip if filename not supplied ---
c
      if( .NOT. is_netcdf_iortem(igrid,1,1) ) goto 9999
c
c   --- set the unit number for surface emissions file ---
c
      iounit = iortem(igrid,1,1)
      fname = temfil(igrid,1,1)
      buffer_offset = buffer_offset_iortem(igrid,1,1)
      write(action,'(2A,I3)') 'Reading the RTRAC NetCDF ',
     &                     'gridded emissions file. Grid: ',igrid
     &                                    
c
c   ---- get the number of layers in this file ---
c
      this_var = 'NLAYS'
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NLAYS', nlays_in)
      if( ierr .NE. NF_NOERR ) goto 7000
c
c   --- allocate the local array ---
c
      allocate( emsgrd(nox-2*buffer_offset,noy-2*buffer_offset,nlays_in) )
c
c   --- if number of layers is not 1 make sure it matches grid ---
c
      if( nlays_in .NE. 1 .AND. nlays_in .NE. nlay(igrid) ) goto 7001
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = nox-2*buffer_offset
      data_start(2) = 1
      data_count(2) = noy-2*buffer_offset
      data_start(3) = 1
      data_count(3) = nlays_in
c
c  ---- get the index for timestep containing this time ----
c
      this_date = ndate
      this_time = ttime
      this_tstep = ncf_get_tstep(iounit,action,
     &              ndate,ttime,this_time_tflag,this_time_etflag,
     &                                                       le1day,.TRUE.)
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  ---- loop over the list of emissions spcies ---
c
      do ispc=1,ntotsp
        this_var = ptname(ispc)
        ierr = nf_inq_varid(iounit, this_var, this_varid)
        if( ierr .NE. NF_NOERR) cycle
        ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                               data_count,emsgrd)
        if( ierr .NE. NF_NOERR) goto 7002
c
        do j=1+buffer_offset,nrow(igrid)-buffer_offset
          jtmp = j-buffer_offset
          do i=1+buffer_offset,ncol(igrid)-buffer_offset
            itmp = i-buffer_offset
            do k=1,nlays_in
c
c  --- convert to emissions time and put into array ---
c
                 emstmp = emsgrd(itmp,jtmp,k)/(60.*dtems)
                 emisrt(i,j,k,ispc) = emisrt(i,j,k,ispc) + emstmp
            enddo
c
c  --- next cell ---
c
          enddo
        enddo
c
c  --- next RTRAC emssions species ---
c
      enddo
c
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Number of layers does not match: '
      write(iout,'(A,I4)') 'User supplied: ',nlays_in
      write(iout,'(A,I4)') 'Value in file: ',nlay(igrid)
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
