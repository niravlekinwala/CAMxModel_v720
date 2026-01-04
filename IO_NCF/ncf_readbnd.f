C**** NCF_READBND
c
      subroutine ncf_readbnd(bndtim,bnddate,did_update)
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
c        bndtim     model simulation time (HHMM)
c        bnddate    model simulation date (YYJJJ)
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
      integer bnddate
      real    bndtim
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
      integer       this_varid, n3d, n4d, i, j, k, ierr
      integer       this_time_tflag, this_time_etflag
      real          date_time, convfac
c
      real, allocatable, dimension(:,:,:) :: bctmp
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      did_update = .FALSE.
      action = 'Reading lateral boundary conditions file.'
c
c  --- check if current timestep still works ---
c
      date_time = REAL(bnddate)+bndtim/2400.
      if( date_time .GE. bnddate_time_tflag .AND. date_time .LT. 
     &                                  bnddate_time_etflag ) goto 9999
c
c  --- reading the file for more data so echo that ---
c
      did_update = .TRUE.
      write(*,'(a20,$)') 'ncf_readbnd ......'
c
c   --- get the timestep that encompasses this date/time ---
c
      this_tstep = ncf_get_tstep(ibc,action,bnddate,bndtim,
     &          bnddate_time_tflag,bnddate_time_etflag,.FALSE.,.TRUE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol(1)
      data_start(2) = 1
      data_count(2) = nrow(1)
      data_start(3) = 1
      data_count(3) = nlay(1)
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  --- allocate temporary array ---
c
      allocate( bctmp(ncol(1),nrow(1),nlay(1)) )
c
c  ---- loop over model species and skip if not in the IC file ---
c
      do ispc=1,nspec
c
c  --- skip if species not in file ---
c
         if( lbcmap(ispc) .LE. 0 ) then
c
c  --- West boundary ---
c
            do k=1,nlay(1)
               do j=2,nrow(1)-1
                 i = 1
                 n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                 n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
                 convfac = 1.
                 if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                 conc(n4d) = convfac*bdnl(ispc)
c
c  --- East boundary ---
c
                 i = ncol(1)
                 n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                 n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
                    convfac = 1.
                 if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                 conc(n4d) = convfac*bdnl(ispc)
               enddo
c
c  --- South boundary ---
c
               do i=2,ncol(1)-1
                 j = 1
                 n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
                 n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
                 convfac = 1.
                 if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                 conc(n4d) = convfac*bdnl(ispc)
c
c  --- North boundary ---
c
                 j = nrow(1)
                 n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
                 n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
                 convfac = 1.
                 if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                 conc(n4d) = convfac*bdnl(ispc)
               enddo
            enddo
            cycle
         endif
c
c  ---- load this species into local array ---
c
         this_varid = lbcmap(ispc)
         ierr = nf_get_vara_real(ibc,this_varid,data_start,data_count,bctmp)
         if( ierr .NE. NF_NOERR) then
           write(iout,'(//,a)') 'ERROR in NCF_READBND:'
           write(iout,*)'Cannot read boundary conditions data for speices: ',
     &                                   spname(ispc)(:istrln(spname(ispc)))
           write(iout,'(A,I5)') 'NetCDF error code: ',ierr
           call camxerr()
         endif
c
c  ---  Load boundary concentrations; convert gasses from ppm to umol/m3,
c       PM stays at ug/m3 ---
c
         do k=1,nlay(1)
c
c  --- West boundary ---
c
            do j=2,nrow(1)-1
              i = 1
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
              convfac = 1.
              if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              conc(n4d) = MAX(bdnl(ispc),bctmp(i,j,k))*convfac
c
c  --- East boundary ---
c
              i = ncol(1)
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
              convfac = 1.
              if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              conc(n4d) = MAX(bdnl(ispc),bctmp(i,j,k))*convfac
            enddo
c
c  --- South boundary ---
c
            do i=2,ncol(1)-1
              j = 1
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
              convfac = 1.
              if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              conc(n4d) = MAX(bdnl(ispc),bctmp(i,j,k))*convfac
c
c  --- South boundary ---
c
              j = nrow(1)
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(ispc - 1)
              convfac = 1.
              if( ispc .LE. ngas ) convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              conc(n4d) = MAX(bdnl(ispc),bctmp(i,j,k))*convfac
            enddo
         enddo
      enddo
c
c  --- write message to out file ---
c
      write(iout,'(a40,f7.0,i8.5)') 'Read boundary condition file at ',
     &                                                  bndtim,bnddate
c
c  --- dealocate the local array ---
c
      deallocate( bctmp )
      write(*,'(a)') '   Done'
      call flush(6)
c
 9999 continue
      return
      end
