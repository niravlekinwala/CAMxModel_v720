c**** NCF_RD_SA_BCFILE
c
      subroutine ncf_rd_sa_bcfile()
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
c     Reads the source apportionment BC file and loads the values into
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
      integer       data_start(4), data_count(4), this_tstep, ispc, itr
      integer       this_varid, date_time_tflag, date_time_etflag, ierr
      integer       n3d, n4d, igrd, ip, ic, ig, i, j, k
      real          convfac
c
      real,         allocatable,  dimension(:,:,:) ::  bndgrd
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
      data_count(3) = nlay(1)
c
c   --- get the timestep that encompasses this date/time ---
c
      action = 'Reading NetCDF SA boundary conditions file.'
      this_tstep = ncf_get_tstep(iorbc,action,begdate,begtim,
     &                 date_time_tflag,date_time_etflag,.FALSE.,.TRUE.)
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  ---- allocate the temp array ---
c
      allocate( bndgrd(ncol(1),nrow(1),nlay(1)) )
c
c  ---- loop over model species in the BC file ---
c
      do ispc=1,num_iorbc
          this_var = spc_iorbc(ispc)
          do itr=1,ntotsp
             if( spc_iorbc(ispc) .NE. ptname(itr) ) cycle
c
c  --- find this species in the file ---- 
c
             ierr = nf_inq_varid(iorbc, this_var, this_varid)
             if( ierr .NE. NF_NOERR) goto 7000
             ierr = nf_get_vara_real(iorbc,this_varid,data_start,
     &                                               data_count,bndgrd)
             if( ierr .NE. NF_NOERR) goto 7001
c
c  ---- put concs into the global array at correct offset ---
c
             do k = 1,nlay(1)
                do j = 2,nrow(1)-1
                  i = 1
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bndgrd(i,j,k))*convfac
                  i = ncol(1)
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bndgrd(i,j,k))*convfac
                enddo

                do i = 2,ncol(1)-1
                  j = 1
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bndgrd(i,j,k))*convfac

                  j = nrow(1)
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bndgrd(i,j,k))*convfac
                enddo
             enddo
          enddo
c
c  --- next species ---
c
      enddo
c
      deallocate( bndgrd )
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 write(iout,'(//,A)') 'ERROR in NCF_RD_SA_BCFILE:'
      write(iout,'(A)') TRIM(action)
      write(iout,'(2A)') 'Cannot find SA BC species in file: ',TRIM(this_var)
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 write(iout,'(//,A)') 'ERROR in NCF_RD_SA_BCFILE:'
      write(iout,'(A)') TRIM(action)
      write(iout,'(2A)') 'Error reading gridded data from SA BC species ',
     &                                    TRIM(this_var),' in SA BC file.'
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
