c**** NCF_GET_NLAYERS
c
      function ncf_get_nlayers(iounit,fname,action)
      use filunit
      use grid
      implicit none
      integer ncf_get_nlayers
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the global file attributes for the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         iounit I  NCF file ID
c         fname  C  filename
c         action C  string that describes file being read
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
      include 'camx.prm'
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer       iounit
      character*(*) fname
      character*(*) action
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
      integer this_dimid, nlays_in, ierr
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- assume one layer ---
c
      ncf_get_nlayers = 1
c
c  --- get the id for the layer dimension ---
c
      ierr = nf_inq_dimid(iounit, "LAY", this_dimid  )
      if( ierr .NE. NF_NOERR ) goto 7000
c
c  --- get the value for this file ---
c
      ierr = nf_inq_dimlen(iounit,this_dimid,nlays_in)
      if( ierr .NE. NF_NOERR ) goto 7000
c
c  --- set return value to this ---
c
      ncf_get_nlayers = nlays_in
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_GET_NLAYERS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for layer dimension.'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
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
