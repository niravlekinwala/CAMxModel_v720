c**** NCF_WRT_TOPO
c
      subroutine ncf_wrt_topo(action,igrd,iounit,grid_ncol,grid_nrow,grid_topo)
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
c           igrd        I grid number
c           iounit      I NetCDF file ID of file
c           grid_ncol   I number of grid cells in X direction
c           grid_nrow   I number of grid cells in Y direction
c           grid_topo   R height above sea level of cell center
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
      integer       igrd
      integer       iounit
      integer       grid_ncol
      integer       grid_nrow
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
      integer       this_varid, icol, irow, ierr, itmp, jtmp
      integer       ibegcel, iendcel, jbegcel, jendcel, buffer_offset
      real,         allocatable, dimension(:,:) :: rarray_2d
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- set cell offset for nested grids --
c
      if( igrd .EQ. 1 ) then
         allocate( rarray_2d(grid_ncol,grid_nrow) )
         ibegcel = 1
         iendcel = grid_ncol
         jbegcel = 1
         jendcel = grid_nrow
         buffer_offset = 0
      else
         allocate( rarray_2d(grid_ncol-2,grid_nrow-2) )
         ibegcel = 2
         iendcel = grid_ncol-1
         jbegcel = 2
         jendcel = grid_nrow-1
         buffer_offset = 1
      endif
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
 
