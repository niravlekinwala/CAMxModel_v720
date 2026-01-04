c**** NCF_WRT_DIM
c
      subroutine ncf_wrt_dim(action,iounit,igrd,numcols,numrows,nlays,nspcs)
      use ncf_iomod
      use filunit
      use grid
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine writes the dimensions to the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c           action  C  name of file to open
c           iounit  I  NetCDF file ID of file
c           igrd    I  grid number of this grid
c           numcols I number of cols in this file
c           numrows I number of cols in this file
c           nlays   I  number of layers in this file
c           nspcs   I  number of species in the file
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c     06/03/21   --gwilson--    Changed var dimension to be model species
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*(*) action
      integer       iounit
      integer       igrd
      integer       numcols
      integer       numrows
      integer       nlays
      integer       nspcs
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
      integer ierr
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ncf_date_time = 2
      ncf_lay = nlays
      ncf_col = numcols
      ncf_row = numrows
      if( igrd .GT. 1 ) then
        ncf_col = numcols - 2
        ncf_row = numrows - 2
      endif
      ncf_var = nspcs+1
c
      ierr = nf_def_dim(iounit, "TSTEP", NF_UNLIMITED, ncf_tstep_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "DATE-TIME", ncf_date_time, ncf_date_time_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "LAY", ncf_lay, ncf_lay_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "VAR", ncf_var, ncf_var_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "ROW", ncf_row, ncf_row_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_dim(iounit, "COL", ncf_col, ncf_col_dimid )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_DIM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot write dimensions to file.'
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
 
