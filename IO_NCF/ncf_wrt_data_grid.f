c**** NCF_WRT_DATA_GRID
c
      subroutine ncf_wrt_data_grid(action,iounit,igrd,grid_ncol,grid_nrow,
     &               grid_xorig,grid_yorig,grid_deltax,grid_deltay,nlays,
     &                                       grid_lat,grid_lon,grid_topo)
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
c           igrd        I grid index for this grid
c           grid_ncol   I number of grid cells in X direction
c           grid_nrow   I number of grid cells in Y direction
c           grid_xorig  R X coordinate of grid origin      
c           grid_yorig  R Y coordinate of grid origin      
c           grid_dxcell R cell width in X direction
c           grid_dycell R cell width in Y direction
c           nlays       I number of layers in this file
c           grid_lat    R latitude value of cell center
c           grid_lon    R longitude value of cell center
c           grid_topo   R height above sea level of cell center
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c     03/15/21   --gwilson--    Changed TOPO to float from double
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
      integer       igrd
      integer       grid_ncol
      integer       grid_nrow
      real          grid_xorig
      real          grid_yorig
      real          grid_deltax
      real          grid_deltay
      integer       nlays
      real          grid_lat(grid_ncol,grid_nrow)
      real          grid_lon(grid_ncol,grid_nrow)
      real          grid_topo(grid_ncol,grid_nrow)
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
      real          factor
      integer,      allocatable, dimension(:)   :: iarray_1d
      integer       this_varid, icol, irow, ilay, ierr, itmp, jtmp
      integer       ibegcel, iendcel, jbegcel, jendcel, buffer_offset
      real*4,       allocatable, dimension(:,:) :: rarray_2d
      real*8,       allocatable, dimension(:)   :: darray_1d
      real*8,       allocatable, dimension(:,:) :: darray_2d
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- variable for X coordinates ---
c
      factor = 1000.
      if( llatlon ) factor = 1.0
      this_var = "X"
      allocate( darray_1d(ncol(igrd)) )
      do icol=1,ncol(igrd)
        darray_1d(icol) = DBLE((grid_xorig + 
     &                           (REAL(icol)-0.5)*grid_deltax)/factor)
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
      allocate( darray_1d(nrow(igrd)) )
      do irow=1,nrow(igrd)
        darray_1d(irow) = DBLE((grid_yorig + 
     &                           (REAL(irow)-0.5)*grid_deltay)/factor)
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
      allocate( iarray_1d(nlays) )
      do ilay=1,nlays
        iarray_1d(ilay) = ilay
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_int(iounit,this_varid,iarray_1d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( iarray_1d )
c
c  --- variable for longitude coordinates ---
c
      if( igrd .EQ. 1 ) then
         allocate( darray_2d(ncol(igrd),nrow(igrd)) )
         allocate( rarray_2d(ncol(igrd),nrow(igrd)) )
         ibegcel = 1
         iendcel = ncol(igrd)
         jbegcel = 1
         jendcel = nrow(igrd)
         buffer_offset = 0
      else
         allocate( darray_2d(ncol(igrd)-2,nrow(igrd)-2) )
         allocate( rarray_2d(ncol(igrd)-2,nrow(igrd)-2) )
         ibegcel = 2
         iendcel = ncol(igrd)-1
         jbegcel = 2
         jendcel = nrow(igrd)-1
         buffer_offset = 1
      endif
c
      this_var = "longitude"
      do icol=ibegcel,iendcel
        itmp = icol-buffer_offset
        do irow=jbegcel,jendcel
           jtmp = irow-buffer_offset
           darray_2d(itmp,jtmp) = DBLE(grid_lon(icol,irow))
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
      do icol=ibegcel,iendcel
        itmp = icol-buffer_offset
        do irow=jbegcel,jendcel
           jtmp = irow-buffer_offset
           darray_2d(itmp,jtmp) = DBLE(grid_lat(icol,irow))
        enddo
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_double(iounit,this_varid,darray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
c
c  --- variable for topo ---
c
      this_var = "topo"
      do icol=ibegcel,iendcel
        itmp = icol-buffer_offset
        do irow=jbegcel,jendcel
           jtmp = irow-buffer_offset
           rarray_2d(itmp,jtmp) = grid_topo(icol,irow)
        enddo
      enddo
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_var_real(iounit,this_varid,rarray_2d)
      if( ierr .NE. NF_NOERR ) goto 7001
      deallocate( rarray_2d )
      deallocate( darray_2d )
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
 
