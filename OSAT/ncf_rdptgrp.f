C***** NCF_RDPTGRP
c
      subroutine ncf_rdptgrp(igroup,numcls,emscls,izcel)
      use chmstry
      use filunit
      use ptemiss
      use tracer
      implicit none
c
c      Copyright 1996 - 2022
c     Ramboll
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
      include 'camx.inc'
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igroup
      integer numcls
      real    emscls(MXTRCLS,*)
      integer izcel(*)
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
      character*200 action, fname
      character*10  this_var
      integer       this_tstep, ierr, ispc, idxfile, iounit
      integer       this_time_tflag, this_time_etflag, num_ptsfiles
      integer       idxpt, this_varid, data_start(2), data_count(2)
      integer       icls, idx_start
c
      real,    allocatable, dimension(:) ::  emispts
      real,    allocatable, dimension(:) ::  emispts_in
      integer, allocatable, dimension(:) ::  sa_region
      real,    allocatable, dimension(:) ::  plume_bot
      real,    allocatable, dimension(:) ::  plume_top
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- skip if filename not supplied ---
c
      num_ptsfiles = num_iortpt(igroup)
      if( igroup .EQ. 0 ) num_ptsfiles = npoint_files
      do idxfile=1,num_ptsfiles
         if( .NOT. ltptfl(igroup,idxfile) ) cycle
         write(action,'(2A,I3,A,I3)') 'Reading the SA NetCDF point source',
     &            ' emissions file for Group: ',igroup,' File: ',idxfile
         idx_start = idx_start_sapts(igroup,idxfile)
c
c   --- set the unit number for file ---
c
         if( igroup .EQ. 0 ) then
             iounit = iptem(idxfile)
             write(fname,'(A,I3)') 'PTSOURCE -- UNIT ',iptem(idxfile)
             if( .NOT. is_netcdf_iptem(idxfile) ) cycle
         else
             if( .NOT. is_netcdf_iortpt(igroup,idxfile) ) cycle
             iounit = iortpt(igroup,idxfile)
             fname = tptfil(igroup,idxfile)
         endif
c
c  ---- allocate the temporary array ---
c
         allocate( emispts(nptsrc_safile(igroup,idxfile)) )
         allocate( emispts_in(nptsrc_safile(igroup,idxfile)) )
         allocate( sa_region(nptsrc_safile(igroup,idxfile)) )
         allocate( plume_bot(nptsrc_safile(igroup,idxfile)) )
         allocate( plume_top(nptsrc_safile(igroup,idxfile)) )
         plume_bot = 0.
         plume_top = 0.
c
c  ---- get the region override flag ---
c
         sa_region = 0
         this_var = 'saoverride'
         ierr = nf_inq_varid(iounit, this_var, this_varid)
         if( ierr .NE. NF_NOERR ) cycle
         ierr = nf_get_var_int(iounit,this_varid,sa_region)
         if( ierr .NE. NF_NOERR ) goto 7001
         do idxpt=1,nptsrc_safile(igroup,idxfile)
            if( sa_region(idxpt) .GT. 0 ) izcel(idxpt+idx_start) = -sa_region(idxpt)
         enddo
c
c  ---- get the index for timestep containing this time ----
c
         this_tstep = ncf_get_tstep(iounit,action,date,time,
     &               this_time_tflag,this_time_etflag,le1day,.TRUE.)
c
c   --- set the indexes for what to read ---
c
         data_start(1) = 1
         data_count(1) = nptsrc_safile(igroup,idxfile)
         data_start(2) = this_tstep
         data_count(2) = 1

         do idxpt = 1,nptsrc_safile(igroup,idxfile)
           effph(idxpt+idx_start) = 0.
           flowrat(idxpt+idx_start) = 0.
         enddo
c
c  ---- plumerise ---
c
         this_var = 'plumerise'
         ierr = nf_inq_varid(iounit, this_var, this_varid)
         if( ierr .EQ. NF_NOERR) then
            ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                  data_count,emispts_in)
            if( ierr .NE. NF_NOERR) goto 7001
c
c  --- Put the data into the global array ----
c
            do idxpt=1,nptsrc_safile(igroup,idxfile)
               if( emispts_in(idxpt) .GT. 0 ) 
     &               effph(idxpt+idx_start) = emispts_in(idxpt)
            enddo
         endif
c
c   --- read the emissions for this hour ---
c
         do ispc=1,nspec
            if( lemmap(ispc) .LE. 0 ) cycle
c
c  --- load data into temporary arrray ---
c
             this_var = emspcname(lemmap(ispc))
             ierr = nf_inq_varid(iounit, this_var, this_varid)
             if( ierr .NE. NF_NOERR) cycle
             ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                               data_count,emispts_in)
             if( ierr .NE. NF_NOERR) goto 7001
c
c  --- put into regular model array ---
c
              do idxpt=1,nptsrc_safile(igroup,idxfile)
                ptemis(idxpt+idx_start,lemmap(ispc)) = emispts_in(idxpt)/(60.*dtems)
              enddo
c
c  --- if not needed for source apportionment skip the rest ----
c
             if( .NOT. lusespc(ispc) ) cycle
c
c   --- put into regular model array ---
c
              do idxpt=1,nptsrc_safile(igroup,idxfile)
                emispts(idxpt) = emispts_in(idxpt)/(60.*dtems)
c
c   --- load into the tracer emissions array ---
c
                do icls=1,ntrcls
c
c   --- if this is the NOx emissions tracer and this is a 
c       PiG source, skip it (PiG treated elsewhere) ---
c
                   if( lpigsa(idxpt+idx_start) ) cycle
                   emscls(icls,idxpt+idx_start) = 
     &                   emscls(icls,idxpt+idx_start) + 
     &                              emispts(idxpt) * trspmap(ispc,icls)
                enddo
              enddo
c
c   --- next species ---
c
         enddo
c
c   --- get plume bottom, if available ---
c
         this_var = 'plume_bottom'
         ierr = nf_inq_varid(iounit, this_var, this_varid)
         if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                        data_count,plume_bot)
             if( ierr .NE. NF_NOERR) goto 7001
         endif
c
c   --- get plume top, if available ---
c
         this_var = 'plume_top'
         ierr = nf_inq_varid(iounit, this_var, this_varid)
         if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                        data_count,plume_top)
             if( ierr .NE. NF_NOERR) goto 7001
         endif
c
c  --- Put the data into the global array ----
c
         do idxpt=1,nptsrc_safile(igroup,idxfile)
            if( plume_top(idxpt) .GT. 0. ) then
               if( plume_bot(idxpt) .GE. plume_top(idxpt) ) goto 7002
               if( plume_bot(idxpt) .LT. 0. ) goto 7002
               effph(idxpt+idx_start) = -plume_bot(idxpt)
               flowrat(idxpt+idx_start) = -plume_top(idxpt)
            endif
         enddo
c
c   --- next file ---
c
         deallocate( emispts    )
         deallocate( emispts_in )
         deallocate( sa_region  )
         deallocate( plume_bot  )
         deallocate( plume_top  )
      enddo
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTGRP: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTGRP: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTGRP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Point source found with invalid plume ',
     &                                        'distribution values.'
      write(iout,'(A,I10)') 'Point source number in file: ',idxpt
      write(iout,'(A,F10.4)') 'Plume bottom: ',plume_bot(idxpt)
      write(iout,'(A,F10.4)') 'Plume top   : ',plume_top(idxpt)
      call camxerr()
c
 9999 continue
c
      return
      end
