c**** NCF_WRT_DATA_GRIDSMP
c
      subroutine ncf_wrt_data_gridsmp(action,iounit,grid_ncol,grid_nrow,
     &                       grid_xorig,grid_yorig,grid_deltax,grid_deltay)
      use ncf_iomod
      use filunit
      use grid
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine writes the variable data for the grid definiton 
c    variables to the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll 
c      Argument description:
c       Inputs:
c           action      C name of file to open
c           iounit      I NetCDF file ID of file
c           grid_ncol   I number of grid cells in X direction
c           grid_nrow   I number of grid cells in Y direction
c           grid_xorig  R X coordinate of grid origin      
c           grid_yorig  R Y coordinate of grid origin      
c           grid_deltax R cell width in X direction
c           grid_deltay R cell width in Y direction
c       Outputs:
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
      include 'netcdf.inc'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*(*) action
      integer       iounit
      integer       grid_ncol
      integer       grid_nrow
      real          grid_xorig
      real          grid_yorig
      real          grid_deltax
      real          grid_deltay
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 this_var
      real*8,       allocatable, dimension(:)   :: darray_1d
      real*8,       allocatable, dimension(:,:) :: darray_2d
      real,         allocatable, dimension(:,:) :: grid_lat
      real,         allocatable, dimension(:,:) :: grid_lon
      real,         allocatable, dimension(:,:) :: grid_mapscl
      real,         allocatable, dimension(:)   :: grid_xsize
      real          grid_ysize
      integer,      allocatable, dimension(:)   :: iarray_1d
      integer       this_varid, icol, irow, ilay, ierr, nlays
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- allocate the arrays for lat/lon for this grid ---
c
      allocate( grid_lat(grid_ncol,grid_nrow) )
      allocate( grid_lon(grid_ncol,grid_nrow) )
      allocate( grid_mapscl(grid_ncol,grid_nrow) )
      allocate( grid_xsize(grid_nrow) )
c
c  --- call routine to fill values ---
c
      if( llatlon ) then
          call grdprep(grid_ncol,grid_nrow,grid_lon,grid_lat,
     &       grid_mapscl,grid_xorig,grid_yorig,grid_deltax,
     &                         grid_deltay,grid_xsize,grid_ysize)
      else
          call grdprep(grid_ncol,grid_nrow,grid_lon,grid_lat,
     &       grid_mapscl,grid_xorig/1000.,grid_yorig/1000.,
     &            grid_deltax/1000.,grid_deltay/1000.,grid_xsize,grid_ysize)
      endif
c
c  --- variable for X coordinates ---
c
      this_var = "X"
      allocate( darray_1d(grid_ncol) )
      do icol=1,grid_ncol
        darray_1d(icol) = DBLE((grid_xorig + 
     &                           (REAL(icol)-0.5)*grid_deltax)/1000.)
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_1d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( darray_1d )
c
c  --- variable for Y coordinates ---
c
      this_var = "Y"
      allocate( darray_1d(grid_nrow) )
      do irow=1,grid_nrow
        darray_1d(irow) = DBLE((grid_yorig + 
     &                           (REAL(irow)-0.5)*grid_deltay)/1000.)
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_1d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( darray_1d )
c
c  --- variable for Layer heights ---
c
      this_var = "layer"
      allocate( iarray_1d(1) )
      iarray_1d(1) = 1
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_int(iounit,this_varid,iarray_1d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( iarray_1d )
c
c  --- variable for longitude coordinates ---
c
      allocate( darray_2d(grid_ncol,grid_nrow) )
c
      this_var = "longitude"
      do icol=1,grid_ncol
        do irow=1,grid_nrow
           darray_2d(icol,irow) = DBLE(grid_lon(icol,irow))
        enddo
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
c
c  --- variable for latitude coordinates ---
c
      this_var = "latitude"
      do icol=1,grid_ncol
        do irow=1,grid_nrow
           darray_2d(icol,irow) = DBLE(grid_lat(icol,irow))
        enddo
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
c
c  --- deallocate local arrays ---
c
      deallocate( darray_2d )
      deallocate( grid_lat )
      deallocate( grid_lon )
      deallocate( grid_mapscl )
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_DATA_GRID:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_DATA_GRID:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot write data for the variable: ',
     &                                      this_var(:istrln(this_var))
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
 
