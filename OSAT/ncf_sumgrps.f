      subroutine ncf_sumgrps(numcols,numrows,nspmod,nsptrac,ndlast,ttlast,
     &                      emstot,emslft,emsbas,emsoth,emssum,lemit)
      use filunit
      use grid
      use chmstry
      use ptemiss
      use tracer
      implicit none
c
c      Copyright 1996 - 2022
c     Ramboll
c
c
c----CAMx v7.20 220430
c
c     
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c    01/15/21  --gwilson-- Fixed bug in reading/accessing 3D files
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
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
      integer ncf_get_tstep
c
c-----------------------------------------------------------------------
c    Argument Declarations:
c-----------------------------------------------------------------------
c
      integer numcols
      integer numrows
      integer nspmod
      integer nsptrac
      integer ndlast
      real    ttlast
      real*8  emstot(numcols,numrows,nspmod)
      real*8  emslft(numcols,numrows,nspmod)
      real    emsbas(nspmod,nsptrac)
      real    emsoth(nspmod,nsptrac)
      real    emssum(nspmod,nsptrac)
      logical lemit(*)
c
c-----------------------------------------------------------------------
c    Local Variables:
c-----------------------------------------------------------------------
c
      character*200 fname, action
      character*10  this_var
      integer       this_date, iounit, igroup, this_tstep
      integer       iseg, idum, ispc, i, j, k, itmp, jtmp
      integer       igrid, idxfile, ierr, nlays_in, this_varid
      integer       this_time_tflag, this_time_etflag, buffer_offset
      integer       this_dimid, numpts, num_emsfiles, num_ptsfiles
      integer       data_start2(2), data_count2(2)
      integer       data_start4(4), data_count4(4)
      real          this_time
      logical       lpass, in_simulation
c
      real,    allocatable, dimension(:,:,:) ::  emstmp
      real,    allocatable, dimension(:,:,:) ::  emsgrd
      real,    allocatable, dimension(:)     ::  emspnt
      integer, allocatable, dimension(:)     ::  sa_region
      integer, allocatable, dimension(:)     ::  idcompact
      integer, allocatable, dimension(:)     ::  izcel
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- loop over all of the groups ----
c
      do 20 igroup=0,ngroup
c
c   --- only process if filename is supplied for group ---
c
         do 25 igrid=1,ngrid
            num_emsfiles = num_iortem(igrid,igroup)
            if( igroup .EQ. 0 ) num_emsfiles = nemiss_files(igrid)
            do 26 idxfile=1,num_emsfiles
              if( .NOT. larsrc) goto 26
c
c   --- set the unit number for surface emissions file ---
c
              lpass = .FALSE.
              if( igroup .EQ. 0 ) then
                  if( .NOT. is_netcdf_iarem(igrid,idxfile) ) goto 26
                  iounit = iarem(igrid,idxfile)
                  write(fname,'(A,I8)') 'EMISSIONS -- UNIT ',
     &                                                iarem(igrid,idxfile)
                  write(action,'(2A,I3,A,I3)') 'Reading the NetCDF ',
     &                              'gridded emissions file. Grid: ',
     &                                           igrid,' File: ',idxfile
                  buffer_offset = buffer_offset_iarem(igrid,idxfile)
              else
                  if( .NOT. is_netcdf_iortem(igrid,igroup,idxfile) ) goto 26
                  iounit = iortem(igrid,igroup,idxfile)
                  fname = temfil(igrid,igroup,idxfile)
                  write(action,'(2A,I3,A,I3,A,I3)') 'Reading the SA NetCDF ',
     &                     'gridded emissions file. Grid: ',igrid,
     &                                    ' Group: ',igroup,' File: ',idxfile
                  buffer_offset = buffer_offset_iortem(igrid,igroup,idxfile)
              endif
c
c   ---- get the number of layers in this file ---
c
              this_var = 'NLAYS'
              ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NLAYS', nlays_in)
              if( ierr .NE. NF_NOERR ) goto 7000
c
c  ---- allocate the temporary array ---
c
              allocate( emstmp(ncol(igrid)-2*buffer_offset,nrow(igrid)-2*buffer_offset,nlays_in) )
              allocate( emsgrd(ncol(igrid),nrow(igrid),nlays_in) )
c
c   --- if number of layers is not 1 make sure it matches grid ---
c
              if( nlays_in .NE. 1 .AND. 
     &                            nlays_in .NE. nlay(igrid) ) goto 7001
c
c   --- set the indexes for what to read ---
c
              data_start4(1) = 1
              data_count4(1) = ncol(igrid)-2*buffer_offset
              data_start4(2) = 1
              data_count4(2) = nrow(igrid)-2*buffer_offset
              data_start4(3) = 1
              data_count4(3) = nlays_in
c
c   --- loop over time in this simulation ---
c
              this_date = date
              this_time = time
              in_simulation = .TRUE.
              do while( in_simulation ) 
c
c  ---- get the index for timestep containing this time ----
c
                 this_tstep = ncf_get_tstep(iounit,action,
     &              this_date,this_time,this_time_tflag,this_time_etflag,
     &                                                         le1day,.TRUE.)
                 data_start4(4) = this_tstep
                 data_count4(4) = 1
c
c  ---- loop over the list of emissions spcies ---
c
                 do ispc=1,nspec
                     this_var = spname(ispc)
                     ierr = nf_inq_varid(iounit, this_var, this_varid)
                     if( ierr .NE. NF_NOERR) cycle
                     ierr = nf_get_vara_real(iounit,this_varid,data_start4,
     &                                               data_count4,emstmp)
                     if( ierr .NE. NF_NOERR) goto 7003
c
                     do k=1,nlays_in
                        do j=1+buffer_offset,nrow(igrid)-buffer_offset
                            jtmp = j-buffer_offset
                            do i=1+buffer_offset,ncol(igrid)-buffer_offset
                              itmp = i-buffer_offset
                              emsgrd(i,j,k) = emstmp(itmp,jtmp,k)
                           enddo
                        enddo
                     enddo
                     call sum1grd(numcols,numrows,ncol(igrid),nrow(igrid),nlays_in,
     &                             nspmod,nsptrac,igroup,igrid,ispc,emssum,emsgrd,
     &                                          emsbas,emsoth,emslft,emstot,lemit)
                 enddo
c
c  --- update time and check if still in simulation ---
c
                 this_time = this_time/100.
                 call uptime(this_time,this_date,dtems)
                 if(this_time .GE. 24.) then
                   this_time  = this_time - 24.
                   this_date = this_date + 1
                   if( MOD(this_date,1000) .GT. 365 ) then
                      if( MOD(INT(this_date/1000),4) .EQ. 0 ) then
                         if( MOD(this_date,1000) .EQ. 367 )
     &                           this_date = (INT(this_date/1000)+1)*1000 + 1
                      else
                         this_date = (INT(this_date/1000)+1)*1000 + 1
                      endif
                   endif
                 endif
                 if( this_date .GT. ndlast ) in_simulation = .FALSE.
                 if( this_date .EQ. ndlast .AND. this_time .GE. ttlast ) 
     &                                           in_simulation = .FALSE.
                 this_time = this_time*100.
            enddo   !--- end while
c
            deallocate( emsgrd )
            deallocate( emstmp )
c
   26       continue 
   25    continue 
c
c   --- only process if filename is supplied for group ---
c
         num_ptsfiles = num_iortpt(igroup)
         if( igroup .EQ. 0 ) num_ptsfiles = npoint_files
         do idxfile=1,num_ptsfiles
            if( .NOT. ltptfl(igroup,idxfile) .OR. .NOT. lptsrc ) cycle
c
c   --- set the unit number for elevated points emissions file ---
c
            lpass = .FALSE.
            if( igroup .EQ. 0 ) then
                if( .NOT. is_netcdf_iptem(idxfile) ) cycle
                iounit = iptem(idxfile)
                write(fname,'(A,I8)') 'PTSOURCE -- UNIT ',iptem(idxfile)
            else
                if( .NOT. is_netcdf_iortpt(igroup,idxfile) ) cycle
                iounit = iortpt(igroup,idxfile)
                fname = tptfil(igroup,idxfile)
            endif
c
c --- get the number of points (used as dimension for columns) ---
c
            ierr = nf_inq_dimid(iounit, "COL", this_dimid )
            if( ierr .NE. NF_NOERR ) goto 7002
            ierr = nf_inq_dimlen(iounit,this_dimid,numpts)
            if( ierr .NE. NF_NOERR ) goto 7002
c
c ---- allocate array for emissions ----
c
            allocate( emspnt(numpts) )
            allocate( sa_region(numpts) )
            allocate( idcompact(numpts) )
            allocate( izcel(numpts) )
            izcel = 0
            do i=1,numpts
               idcompact(i) = i
            enddo
c
c  --- get the region override for the sources ---
c
            this_var = "saoverride"
            ierr = nf_inq_varid(iounit,this_var,this_varid)
            if( ierr .NE. NF_NOERR ) goto 7004
            ierr = nf_get_var_int(iounit,this_varid,sa_region)
            if( ierr .NE. NF_NOERR ) goto 7003
            izcel = 0
            do i=1,numpts
              if( sa_region(i) .GT. 0 ) izcel(i) = -sa_region(i)
            enddo
c
c   --- loop over time in this simulation ---
c
            this_date = date
            this_time = time
            in_simulation = .TRUE.
            do while( in_simulation )
c
c  ---- get the index for timestep containing this time ----
c
               this_tstep = ncf_get_tstep(iounit,action,
     &              this_date,this_time,this_time_tflag,this_time_etflag,
     &                                                       le1day,.TRUE.)
               data_start2(1) = 1
               data_count2(1) = numpts
               data_start2(2) = this_tstep
               data_count2(2) = 1
c
c  ---- loop over the list of emissions spcies ---
c
               do ispc=1,nspec
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
                 this_var = spname(ispc)
                 ierr = nf_inq_varid(iounit, this_var, this_varid)
                 if( ierr .NE. NF_NOERR) cycle
                 ierr = nf_get_vara_real(iounit,this_varid,data_start2,
     &                                               data_count2,emspnt)
                 call sum1pnt(numcols,numrows,nspmod,nsptrac,igroup,idxfile,
     &                  ispc,nptsrc,emsbas,emsoth,emslft,emstot,emspnt,
     &                                       emssum,izcel,idcompact,lemit)
               enddo
c
c  --- update time and check if still in simulation ---
c
               call uptime(this_time,this_date,dtems)
               if(this_time .GE. 24.) then
                 this_time  = this_time - 24.
                 this_date = this_date + 1
                 if( MOD(this_date,1000) .GT. 365 ) then
                    if( MOD(INT(this_date/1000),4) .EQ. 0 ) then
                       if( MOD(this_date,1000) .EQ. 367 )
     &                           this_date = (INT(this_date/1000)+1)*1000 + 1
                    else
                       this_date = (INT(this_date/1000)+1)*1000 + 1
                    endif
                 endif
               endif
               if( this_date .GT. ndlast ) in_simulation = .FALSE.
               if( this_date .EQ. ndlast .AND. this_time .GE. ttlast )
     &                                           in_simulation = .FALSE.
            enddo   !--- end while
c
c  --- deallocate local arrays ---
c
            deallocate( emspnt )
            deallocate( sa_region )
            deallocate( idcompact )
            deallocate( izcel )
c
c  --- next file ---
c
         enddo
c
c  --- get the next group ---
c
  20  continue
      return
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_SUMGRPS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_SUMGRPS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Number of layers does not match: '
      write(iout,'(A,I4)') 'User supplied: ',nlays_in
      write(iout,'(A,I4)') 'Value in file: ',nlay(igrid)
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_SUMGRPS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for number ',
     &                                      'of point sources (COL)'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in NCF_SUMGRPS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7004 continue
      write(*,'(//,a)') 'ERROR in NCF_SUMGRPS:'
      write(*,'(A)') action(:istrln(action))
      write(*,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
      end
