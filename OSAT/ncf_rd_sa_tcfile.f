c**** NCF_RD_SA_TCFILE
c
      subroutine ncf_rd_sa_tcfile()
      use grid
      use camxfld
      use filunit
      use tracer
      implicit none
c     
c----CAMx v7.20 220430
c
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     Reads the source apportionment TC file and loads the values into
c     the tracer concentration array. This version is for NetCDF files.
c 
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c       Outputs:
c           
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include  files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.inc'
      include 'chmdat.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer ncf_get_tstep
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      character*10  this_var
      integer       data_start(4), data_count(4), this_tstep, ispc
      integer       this_varid, date_time_tflag, date_time_etflag, ierr
      integer       n3d, i, j, k
      real          convfac
c
      real,         allocatable,  dimension(:,:,:) ::  tinit
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      data_start(1) = 1
      data_count(1) = ncol(1)
      data_start(2) = 1
      data_count(2) = nrow(1)
      data_start(3) = 1
      data_count(3) = 1
c
c   --- get the timestep that encompasses this date/time ---
c
      action = 'Reading NetCDF SA initial conditions file.'
      this_tstep = ncf_get_tstep(iortc,action,begdate,begtim,
     &                 date_time_tflag,date_time_etflag,.FALSE.,.TRUE.)
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  ---- allocate the temp array ---
c
      allocate( tinit(ncol(1),nrow(1),1) )
c
c  ---- loop over model species in the TC file ---
c
      do ispc=1,num_iorbc
          this_var = spc_iorbc(ispc)
c
c  --- find this species in the file ---- 
c
          ierr = nf_inq_varid(iortc, this_var, this_varid)
          if( ierr .NE. NF_NOERR) goto 7000
          ierr = nf_get_vara_real(iortc,this_varid,data_start,
     &                                               data_count,tinit)
          if( ierr .NE. NF_NOERR) goto 7001
c
c  ---- put concs into the global array at correct offset ---
c
          do k=1,ncol(1)
            do j = 1,nrow(1)
               n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(ispc - 1)
               ptconc(n3d) = BNDLPT
               ptconc(n3d) = AMAX1(ptconc(n3d),tinit(i,j,1))
            enddo
          enddo
c
c  --- next species ---
c
      enddo
c
      deallocate( tinit )
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 write(iout,'(//,A)') 'ERROR in NCF_RD_SA_TCFILE:'
      write(iout,'(A)') TRIM(action)
      write(iout,'(2A)') 'Cannot find SA TC species in file: ',TRIM(this_var)
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 write(iout,'(//,A)') 'ERROR in NCF_RD_SA_TCFILE:'
      write(iout,'(A)') TRIM(action)
      write(iout,'(2A)') 'Error reading gridded data from SA TC species ',
     &                                    TRIM(this_var),' in SA TC file.'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
