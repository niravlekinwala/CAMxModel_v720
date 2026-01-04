C*** NCF_LUSEPREP
c
      subroutine ncf_luseprep(igrd)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c     NCF_METPREP reads the landuse NETCDF file and checks
c     the global varibles for consistency with user defined values.
c                          
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         igrd    I grid index
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
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrd
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
      character*200 action
      character*10  namelu, name_in, this_type
      character*10  this_var
      integer       ierr, iunit, itmp
c
c-----------------------------------------------------------------------
c    Data Statements:
c-----------------------------------------------------------------------

      data namelu /'LANDUSE   '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      iunit = isurf(igrd)
      this_type = namelu
c
c --- set the string for error messages ---
c
      write(action,'(A,I5)') 'Reading the landuse file for grid: ',igrd
c
c --- get the type of file to make sure it is the correct file type ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iunit, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iunit, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7002
      endif
      if( name_in(:7) .NE. this_type(:7) ) goto 7001
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(iunit,action,igrd,.FALSE.,.TRUE.,.FALSE.,.TRUE.,itmp)
c
      return
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_METPREP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get necessary global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_METPREP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',this_type(:istrln(this_type))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in NCF_METPREP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)',ERR=9999) 'Cannot find global variable for ',
     &                                       'type of file: CAMx_NAME'
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
