c***** NCF_RDARDDM.F
c
      subroutine ncf_rdarddm(igrid,ndate,ttime,nox,noy,nlay_ems,nddmspc,emisddm)
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
c   This routine reads one hour of emissions for the DDM process
c   and fills the approproate arrays.  The emissions file for one grid
c   but each emissions groups is read. This version is for NeetCDF files.
c    Argument descriptions:
c     Outputs:
c       emisddm   R    array of emissions for DDM speices
c     Inputs:
c       igrid     I    grid number
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       nox       I    number of columns in grid
c       noy       I    number of rows in grid
c       nlay_ems  I    number of layers in emissions files
c       nddmspc   I    number of DDM species
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
      integer nddmspc
      real    emisddm(nox,noy,nlay_ems,nddmspc)
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname, action
      character*10  this_var
      integer      nedge, igroup, idxfile, iounit, ierr, itmp,  jtmp
      integer      nlays_in, ispc, data_start(4), data_count(4), this_date
      integer      this_time_tflag, this_time_etflag, this_tstep
      integer      iddm, i, j, k, iptr, imap,  this_varid, buffer_offset
      real         emstmp, this_time
      logical      luse
c
      real, allocatable, dimension(:,:,:) :: emsgrd
      real, allocatable, dimension(:,:,:) :: tmpgrd
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c  --- if not doing emissions groups, just return here ---
c
      if( nemddm .EQ. 0 ) goto 9999
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
      else
        nedge = 1
      endif
c
c   --- loop over all of the groups ----
c
      do igroup=1,ngroup
        do idxfile=1,num_iortem(igrid,igroup)
c
c   --- skip if filename not supplied ---
c
          if( .NOT. is_netcdf_iortem(igrid,igroup,idxfile) ) cycle
c
c   --- set the unit number for surface emissions file ---
c
          iounit = iortem(igrid,igroup,idxfile)
          fname = temfil(igrid,igroup,idxfile)
          write(action,'(2A,I3,A,I3,A,I3)') 'Reading the DDM NetCDF ',
     &                     'gridded emissions file. Grid: ',igrid,
     &                                    ' Group: ',igroup,' File: ',idxfile
          buffer_offset = buffer_offset_iortem(igrid,igroup,idxfile)
c
c   ---- get the number of layers in this file ---
c
          this_var = 'NLAYS'
          ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NLAYS', nlays_in)
          if( ierr .NE. NF_NOERR ) goto 7000
c
c   --- allocate the local array ---
c
          allocate( tmpgrd(nox-2*buffer_offset,noy-2*buffer_offset,nlays_in) )
          allocate( emsgrd(nox,noy,nlays_in) )
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
          do ispc=1,nspec
            this_var = spname(ispc)
            ierr = nf_inq_varid(iounit, this_var, this_varid)
            if( ierr .NE. NF_NOERR) cycle
            ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                               data_count,tmpgrd)
            if( ierr .NE. NF_NOERR) goto 7002
c
             do k=1,nlays_in
                do j=1+buffer_offset,nrow(igrid)-buffer_offset
                   jtmp = j-buffer_offset
                   do i=1+buffer_offset,ncol(igrid)-buffer_offset
                      itmp = i-buffer_offset
                      emsgrd(i,j,k) = tmpgrd(itmp,jtmp,k)
                   enddo
                enddo
             enddo
c
c  --- find this species in the modeled species list ---
c
             luse = .FALSE.
             do iddm=1,nemddm
                luse = .FALSE.
                if( emddmsp(iddm) .EQ. spname(ispc) ) luse = .TRUE.
                if( emddmsp(iddm) .EQ. NAMALL  ) luse = .TRUE.
                if( emddmsp(iddm) .EQ. NAMVOC
     &                               .AND. lvocsp(ispc) ) luse = .TRUE.
                if( emddmsp(iddm) .EQ. NAMNOX
     &                               .AND. lnoxsp(ispc) ) luse = .TRUE.
                if( emddmsp(iddm) .EQ. NAMHRV
     &                               .AND. lhrvoc(ispc) ) luse = .TRUE.
c
c  --- if this DDM species matches this modeled species, load it ---
c
                if( luse ) then
                   do j=1,noy-1
                     do i=1,nox-1
c
                       imap = igrmap(0,1,igrid,i,j)
                       if( imap .LE. 0 .OR. imap .GT. nregin ) cycle
                       do k=1,nlays_in
c
c  --- convert to emissions time and put into array ---
c
                            emstmp = emsgrd(i,j,k)/(60.*dtems)
                            iptr = iptddm(ispc) + nicddm +
     &                             nbcddm*nedge + (iddm-1)*ngroup*nregin +
     &                                      (    imap-1)*ngroup + igroup - 1
                            emisddm(i,j,k,iptr) = emisddm(i,j,k,iptr) + emstmp
                       enddo
c
c  --- next cell ---
c
                     enddo
                   enddo
                endif
c
c  --- next DDM emssions species ---
c
             enddo
c
c  --- next model species ---
c
          enddo
c
c  --- next file ---
c
          deallocate( emsgrd )
          deallocate( tmpgrd )
        enddo
c
c  --- get the next group ---
c
       enddo
       goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARDDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARDDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Number of layers does not match: '
      write(iout,'(A,I4)') 'User supplied: ',nlays_in
      write(iout,'(A,I4)') 'Value in file: ',nlay(igrid)
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARDDM:'
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
