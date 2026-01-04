C**** NCF_READCNC
c
      subroutine ncf_readcnc()
      use camxfld
      use filunit
      use grid
      use chmstry
      implicit none
c 
c     
c----CAMx v7.20 220430
c
c     This routine reads the NetCDF initial conditions file and 
c     puts the concentrations in the gridded array for instantaneous
c     concentrations.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c  
c     Input arguments: 
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c             
c     Called by: 
c        CAMx
c
      include 'camx.prm'
      include 'camx.inc'
      include 'netcdf.inc'
c
      integer ncf_get_tstep
c
      character*200 action
      integer       data_start(4), data_count(4), ispc, i, j, k, l
      integer       this_varid, iunit, nx, ny, nz, n3d, n4d, igrd
      integer       ierr, ig, ip, ic, this_tstep
      integer       iicdate_time_tflag, iicdate_time_etflag
      real          convfac
c
      real, allocatable, dimension(:,:,:) :: cinit
c
      integer istrln
c
c-----Entry point
c
      iunit = iic
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
      data_start(1) = 1
      data_count(1) = nx
      data_start(2) = 1
      data_count(2) = ny
      data_start(3) = 1
      data_count(3) = nz
c
c   --- get the timestep that encompasses this date/time ---
c
      if( is_netcdf_iic ) then
         action = 'Reading initial conditions file.'
         this_tstep = ncf_get_tstep(iic,action,begdate,begtim,
     &          iicdate_time_tflag,iicdate_time_etflag,.FALSE.,.TRUE.)
         data_start(4) = this_tstep
         data_count(4) = 1
      endif
c
c  ---- allocate the temp array ---
c
      allocate( cinit(nx,ny,nz) )
c
c  ---- loop over model species and skip if not in the IC file ---
c
      do ispc=1,nspec
c
c  --- initialize to lower bound ---
c
         do k=1,nz
           do j = 1,NY
             do i = 1,nx
                n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                n4d = n3d + nx*ny*nz*(ispc - 1)
                conc(n4d) = bdnl(ispc)
             enddo
           enddo
         enddo
c
c  --- skip if species not in file ---
c
         if( licmap(ispc) .LE. 0 ) cycle
         this_varid = licmap(ispc)
c
c  ---- load this species into local array ---
c
         ierr = nf_get_vara_real(iunit,this_varid,data_start,data_count,cinit)
         if( ierr .NE. NF_NOERR) then
           write(iout,'(//,a)') 'ERROR in NCF_READCNC:'
           write(iout,*)'Cannot read initial conditions data for speices: ',
     &                                   spname(ispc)(:istrln(spname(ispc)))
           write(iout,'(A,I5)') 'NetCDF error code: ',ierr
           call camxerr()
         endif
c
c  ---- put concs into the global array at correct offset --- 
c
         do k=1,nz
           do j = 1,ny
             do i = 1,nx
                n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                n4d = n3d + nx*ny*nz*(ispc - 1)
                conc(n4d) = bdnl(ispc)
                conc(n4d) = amax1(conc(n4d),cinit(i,j,k))
             enddo
           enddo
         enddo
      enddo
c
c----- Interpolate coarse grid concentrations to all fine grids ----
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
        nx = ncol(igrd)
        ny = nrow(igrd)
        nz = nlay(igrd)
        do l = 1,nspec
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                n4d = n3d + nx*ny*nz*(l - 1)
                if (l.le.ngas) then
                  convfac = densfac*273./tempk(iptr3d(igrd)-1+n3d)*
     &                        press(iptr3d(igrd)-1+n3d)/1013.
                else
                  convfac = 1.
                endif
                conc(iptr4d(igrd)-1+n4d) = convfac*
     &                     AMAX1( bdnl(l), conc(iptr4d(igrd)-1+n4d) )
              enddo
            enddo
          enddo
        enddo
      enddo
c
c  ---- write message to out file ---
c
      write(iout,'(a40,f7.0,i8.5)') 'Read initial condition file at ',
     &                                                  begtim,begdate
c
c  ---- deallocate the temp array ---
c
      deallocate( cinit )
c
      return
      end
