C***** RDARGRP
c
      subroutine ncf_rdargrp(igrid,ndate,ttime,nox,noy,nlay_ems,
     &                                          numcls,igroup,emscls)
      use filunit
      use grid
      use chmstry
      use bndary
      use camxcom
      use tracer
      implicit none
c
c      Copyright 1996 - 2022
c     Ramboll
c
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      integer ndate
      real    ttime
      integer nox
      integer noy
      integer nlay_ems
      integer numcls
      integer igroup
      real    emscls(numcls,nox,noy,nlay_ems)
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
      character*200 fname, action
      character*10  this_var
      character*4   iname(10)
      integer       ispc, num_emsfiles, idxfile, this_dat
      integer       iounit, iseg, i, j, k, nlays_in, buffer_offset
      integer       data_start(4), data_count(4), itmp, jtmp
      integer       icls, ierr, this_date, this_varid
      integer       this_time_tflag, this_time_etflag, this_tstep
      real          this_time
      logical       lfound, use_this_file
c
      real, allocatable, dimension(:,:,:) :: emsgrd
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- skip if filename not supplied ---
c
      if( .NOT. larsrc ) goto 9999
      num_emsfiles = num_iortem(igrid,igroup)
      if( igroup .EQ. 0 ) num_emsfiles = nemiss_files(igrid)
c
c   --- loop over all files ---
c
      do 10 idxfile=1,num_emsfiles
        use_this_file = .TRUE.
c
c   --- set the unit number for surface emissions file ---
c
        if( igroup .EQ. 0 ) then
            iounit = iarem(igrid,idxfile)
            write(fname,'(A,I10)') 'EMISSIONS -- UNIT ',iarem(igrid,idxfile)
            if( .NOT. is_netcdf_iarem(igrid,idxfile) ) use_this_file = .FALSE.
            write(action,'(2A,I3,A,I3)') 'Reading the NetCDF ',
     &                              'gridded emissions file. Grid: ',
     &                                           igrid,' File: ',idxfile
            buffer_offset = buffer_offset_iarem(igrid,idxfile)
        else
            iounit = iortem(igrid,igroup,idxfile)
            fname = temfil(igrid,igroup,idxfile)
            if( .NOT. is_netcdf_iortem(igrid,igroup,idxfile) ) use_this_file = .FALSE.
            write(action,'(2A,I3,A,I3,A,I3)') 'Reading the SA NetCDF ',
     &                     'gridded emissions file. Grid: ',igrid,
     &                                    ' Group: ',igroup,' File: ',idxfile
            buffer_offset = buffer_offset_iortem(igrid,igroup,idxfile)
        endif
        if( .NOT. use_this_file ) cycle
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
        do ispc=1,nspec
            this_var = spname(ispc)
            ierr = nf_inq_varid(iounit, this_var, this_varid)
            if( ierr .NE. NF_NOERR) cycle
            ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                               data_count,emsgrd)
             if( ierr .NE. NF_NOERR) goto 7002
c
c   --- skip if the species is a not used for SA ---
c
            if( .NOT. lusespc(ispc) ) cycle
c
c   --- convert to correct units ---
c
            do j=1,noy-2*buffer_offset
              do i=1,nox-2*buffer_offset
                do k=1,nlays_in
                   emsgrd(i,j,k) = emsgrd(i,j,k)/(60.*dtems)
                enddo
              enddo
            enddo
c
c   --- load the species into the cells in tracer array ----
c
             do k=1,nlays_in
                do j=1+buffer_offset,noy-buffer_offset
                   jtmp = j-buffer_offset
                   do i=1+buffer_offset,nox-buffer_offset
                     itmp = i-buffer_offset
                     do icls=1,ntrcls
                           emscls(icls,i,j,k) = emscls(icls,i,j,k) + 
     &                              emsgrd(itmp,jtmp,k) * trspmap(ispc,icls)
                    enddo
                  enddo
                enddo
             enddo
c
c   --- next species ---
c
        enddo
c
c  --- next file in grid/group ---
c
      deallocate( emsgrd )
 10   continue
c
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARGRP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARGRP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Number of layers does not match: '
      write(iout,'(A,I4)') 'User supplied: ',nlays_in
      write(iout,'(A,I4)') 'Value in file: ',nlay(igrid)
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDARGRP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
c
      return
      end
