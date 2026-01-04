c**** NCF_RD_SA_ICFILE
c
      subroutine ncf_rd_sa_icfile()
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
c     Reads the source apportionment IC file and loads the values into
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
      real,         allocatable,  dimension(:,:,:) ::  cinit
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
      action = 'Reading NetCDF SA initial conditions file.'
      this_tstep = ncf_get_tstep(ioric,action,begdate,begtim,
     &                 date_time_tflag,date_time_etflag,.FALSE.,.TRUE.)
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  ---- allocate the temp array ---
c
      allocate( cinit(ncol(1),nrow(1),nlay(1)) )
c
c  ---- loop over model species in the IC file ---
c
      do ispc=1,num_ioric
          this_var = spc_ioric(ispc)
          do itr=1,ntotsp
            if( spc_ioric(ispc) .NE. ptname(itr) ) cycle
c
c  --- find this species in the file ---- 
c
             ierr = nf_inq_varid(ioric, this_var, this_varid)
             if( ierr .NE. NF_NOERR) goto 7000
             ierr = nf_get_vara_real(ioric,this_varid,data_start,
     &                                               data_count,cinit)
             if( ierr .NE. NF_NOERR) goto 7001
c
c  ---- put concs into the global array at correct offset ---
c
             do k=1,nlay(1)
               do j = 2,nrow(1)-1
                 do i = 2,ncol(1)-1
                    n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                    n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                    ptconc(n4d) = BNDLPT
                    ptconc(n4d) = AMAX1(ptconc(n4d),cinit(i,j,k))
                 enddo
               enddo
             enddo
          enddo
c
c  --- next species ---
c
      enddo
c
c  --- Interpolate coarse grid concentrations to all fine grids ----
c
      if( ngrid .GT. 1 ) then
          do ip = 1,ngrid
            do ic = 1,nchdrn(ip)
              ig = idchdrn(ic,ip)
              call intrpcnc(nspec,ncol(ip),nrow(ip),nlay(ip),i1(ig),
     &                      j1(ig),nmesh(ig),ncol(ig),nrow(ig),nlay(ig),
     &                      conc(iptr4d(ip)),conc(iptr4d(ig)) )
            enddo
          enddo
      endif
c
c-----Convert from ppm to umol/m3
c
      do igrd = 1,ngrid
        do ispc = 1,num_ioric
          do itr=1,ntotsp
            if( spc_ioric(ispc) .NE. ptname(itr) ) cycle
            if( .NOT. lsagas(itr) ) cycle
            do k = 1,nlay(igrd)
              do j = 1,nrow(igrd)
                do i = 1,ncol(igrd)
                  n3d = i + ncol(igrd)*(j - 1) + ncol(igrd)*nrow(igrd)*(k - 1)
                  n4d = n3d + ncol(igrd)*nrow(igrd)*nlay(igrd)*(itr - 1)
                  convfac = densfac*273./tempk(iptr3d(igrd)-1+n3d)*
     &                          press(iptr3d(igrd)-1+n3d)/1013.
                  ptconc(ipsa3d(igrd)-1+n4d) = convfac*
     &                       AMAX1(BNDLPT, ptconc(ipsa3d(igrd)-1+n4d) )
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      deallocate( cinit )
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 write(iout,'(//,A)') 'ERROR in NCF_RD_SA_ICFILE:'
      write(iout,'(A)') TRIM(action)
      write(iout,'(2A)') 'Cannot find SA IC species in file: ',TRIM(this_var)
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 write(iout,'(//,A)') 'ERROR in NCF_RD_SA_ICFILE:'
      write(iout,'(A)') TRIM(action)
      write(iout,'(2A)') 'Error reading gridded data from SA IC species ',
     &                                    TRIM(this_var),' in SA IC file.'
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
