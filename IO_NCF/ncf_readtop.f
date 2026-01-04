C**** NCF_READTOP
c
      subroutine ncf_readtop(toptim,topdate,did_update)
      use chmstry
      use filunit
      use camxfld
      use grid
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Arguments: 
c        toptim            model simulation time (HHMM)
c        topdate           model simulation date (YYJJJ)
c        did_update .TRUE. if new data was read
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
      include 'camx.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer topdate
      real    toptim
      logical did_update
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
      integer       this_tstep, data_start(4), data_count(4), ispc
      integer       this_varid, n3d, i, j, ierr, ic, ip, ig
      real          date_time
c
      real, allocatable, dimension(:,:,:) :: tctmp
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      did_update = .FALSE.
      action = 'Reading top boundary conditions file.'
c
c  --- check if current timestep still works ---
c
      date_time = REAL(topdate)+toptim/2400.
      if( date_time .GE. topdate_time_tflag .AND. date_time 
     &                             .LT. topdate_time_etflag ) goto 9999
c
c  --- reading the file for more data so echo that ---
c
      did_update = .TRUE.
      write(*,'(a20,$)') 'ncf_readtop ......'
c
c   --- get the timestep that encompasses this date/time ---
c
      this_tstep = ncf_get_tstep(itc,action,topdate,toptim,
     &        topdate_time_tflag,topdate_time_etflag,.FALSE.,.TRUE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol(1)
      data_start(2) = 1
      data_count(2) = nrow(1)
      data_start(3) = 1
      data_count(3) = 1
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  --- allocate temporary array ---
c
      allocate( tctmp(ncol(1),nrow(1),1) )
c
c  ---- loop over model species and skip if not in the IC file ---
c
      do ispc=1,nspec
c
c  ---  if species not in file, set to lower bound ---
c
         if( ltcmap(ispc) .LE. 0 ) then
c
c --- Load top concentrations for master grid ---
c
            do j=1,nrow(1)
              do i=1,ncol(1)
                 n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(ispc - 1)
                 ctop(n3d) = bdnl(ispc)
              enddo
            enddo
            cycle
         endif
c
c  ---- load this species into local array ---
c
         this_varid = ltcmap(ispc)
         ierr = nf_get_vara_real(itc,this_varid,data_start,data_count,tctmp)
         if( ierr .NE. NF_NOERR) then
           write(iout,'(//,a)') 'ERROR in NCF_READTOP:'
           write(iout,*)'Cannot read boundary conditions data for speices: ',
     &                                   spname(ispc)(:istrln(spname(ispc)))
           write(iout,'(A,I5)') 'NetCDF error code: ',ierr
           call camxerr()
         endif
c
c --- Load top concentrations for master grid ---
c
         do j=1,nrow(1)
           do i=1,ncol(1)
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(ispc - 1)
              ctop(n3d) = bdnl(ispc)
              if( tctmp(i,j,1) .GT. bdnl(ispc) ) ctop(n3d) = tctmp(i,j,1)
            enddo
         enddo
      enddo
c
c  --- dealocate the local array ---
c
      deallocate( tctmp )
c
c --- interpolate topcon to all nested grids
c
      if( ngrid .GT. 1) then
        do ip = 1,ngrid
          do ic = 1,nchdrn(ip)
            ig = idchdrn(ic,ip)
            call intrpcnc(nspec,ncol(ip),nrow(ip),1,i1(ig),
     &                    j1(ig),nmesh(ig),ncol(ig),nrow(ig),1,
     &                    ctop(iptr1lay(ip)),ctop(iptr1lay(ig)))
          enddo
        enddo
      endif
c
c  ---- write message to out file ---
c
      write(iout,'(a40,f7.0,i8.5)') 'Read boundary condition file at ',
     &                                                   toptim,topdate
c
      write(*,'(a)') '   Done'
      call flush(6)
c
 9999 continue
      return
      end
