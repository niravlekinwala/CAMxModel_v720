c**** NCF_AREAPREP
c
      subroutine ncf_areaprep(action,iunit,igrd,
     &                     begtim,begdate,endtim,enddate,buffer_offset)
      use chmstry
      use bndary
      use grid
      use filunit
      implicit none
c 
c----CAMx v7.20 220430
c 
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_AREAPREP reads the NETCDF gridded emissions files and checks
c     the global variables for consistency with user defined values.
c     It also creates the global species list.
c                           
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         action        C string that describes what we are doing
c         iunit         I NetCDF file ID
c         igrd          I grid index
c         begtim        R model begin time
c         begdat        I model begin date (YYJJJ)
c         endtim        R model end time
c         enddate       I model end date (YYJJJ)
c         buffer_offset I cell offset for buffer cells in file (0 or 1)
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
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*200 action
      integer       iunit
      integer       igrd
      real          begtim
      integer       begdate
      real          endtim
      integer       enddate
      integer       buffer_offset
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
      character*10  emfil, name_in, this_var
      integer       ierr, i, j, n
      logical       lno_buffer_cells
c
c-----------------------------------------------------------------------
c    Data Statements:
c-----------------------------------------------------------------------
c
      data emfil /'EMISSIONS '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- get the type of file to make sure it is a emissions file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iunit, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iunit, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7002
      endif
      if( name_in(:9) .NE. emfil(:9) ) goto 7001
c
c --- call routine to make sure grid defintion is consistent ---
c
      lno_buffer_cells = .TRUE.
      if( igrd .EQ. 1 ) lno_buffer_cells = .FALSE.
      call ncf_chk_griddef(iunit,action,igrd,.FALSE.,.TRUE.,
     &                         lno_buffer_cells,.FALSE.,buffer_offset)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iunit,action,begdate,begtim,enddate,endtim,le1day)
c
c --- call routine to setup species mappping array ---
c
      call ncf_set_emiss_mapping(iunit,action)
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_AREAPREP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get necessary global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_AREAPREP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',emfil(:istrln(emfil))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in NCF_AREAPREP:'
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
