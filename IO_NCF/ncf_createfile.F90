!**** NCF_CREATEFILE
!
      subroutine ncf_createfile(fname,action,iounit)
      use camx_includes
      use ncf_iomod
!
!----CAMx v7.20 220430
!
!-----------------------------------------------------------------------
!    Description:
!-----------------------------------------------------------------------
!
!   This routine opens a NetCDF file.
!
!      Copyright 1996 - 2022
!     Ramboll
!      Argument description:
!       Inputs:
!           fname  C  name of file to open
!           action C  description of exactly what this call is doing
!       Outputs:
!           iounit I  NetCDF file ID of opened file
!
!-----------------------------------------------------------------------
!    LOG:
!-----------------------------------------------------------------------
!
!     02/20/17   --gwilson--    Original development
!
!-----------------------------------------------------------------------
!    Include files:
!-----------------------------------------------------------------------
!
      include 'netcdf.inc'
!
!-----------------------------------------------------------------------
!    Argument declarations:
!-----------------------------------------------------------------------
!
      character*(*) fname
      character*(*) action
      integer       iounit
!
!-----------------------------------------------------------------------
!    External functions:
!-----------------------------------------------------------------------
!
      integer istrln
!
!-----------------------------------------------------------------------
!    Local variables:
!-----------------------------------------------------------------------
!
      integer ierr, ncf_format
!
!-----------------------------------------------------------------------
!    Entry point:
!-----------------------------------------------------------------------
!
      ncf_format = NF_FORMAT_CLASSIC
      if( ncf_compress ) ncf_format = OR(NF_NETCDF4,NF_FORMAT_CLASSIC)
      ierr = nf_create(fname, ncf_format, iounit)      
      if( ierr .NE. nf_noerr ) then
         write(*,'(//,a)') 'ERROR in NCF_CREATEFILE:'
         write(*,'(A)') action(:istrln(action))
         write(*,'(2A)') 'Could not open file: ', &
                                         fname(:istrln(fname))
         write(*,'(A,I5)') 'NetCDF error code: ',ierr
         call camxerr()
      endif
!
#ifdef CHUNK
      if( ncf_compress ) then
         ierr = nf_set_chunk_cache(NCF_CACHESIZE, NCF_NELEMS, NCF_PREEMPTION)
         if( ierr .NE. nf_noerr ) then
            write(*,'(//,a)') 'ERROR in NCF_CREATEFILE:'
            write(*,'(A)') action(:istrln(action))
            write(*,'(2A)') 'Could not set chunk parameters for: ', &
                                          fname(:istrln(fname))
            write(*,'(A,I5)') 'NetCDF error code: ',ierr
            call camxerr()
         endif
      endif
#endif
!
      goto 9999
!
!-----------------------------------------------------------------------
!    Return point:
!-----------------------------------------------------------------------
!
 9999 continue
      return
      end
