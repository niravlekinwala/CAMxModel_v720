c**** NCF_TOPPREP
c
      subroutine ncf_topprep(begtim,begdate,endtim,enddate)
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
c     NCF_TOPPREP reads the boundary conditions NETCDF file and checks
c     the global varibles for consistency with user defined values.
c     It also prepares mapping variables that map the TC species list
c     to the internal CAMx species list
c                           
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         begtim  R model begin time
c         begdat  I model begin date (YYJJJ)
c         endtim  R model end time
c         enddate I model end date (YYJJJ)
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
      real    begtim
      integer begdate
      real    endtim
      integer enddate
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
      character*10  tcfil, name_in, this_var
      integer       ierr, i, j, n, itmp
c
c-----------------------------------------------------------------------
c    Data Statements:
c-----------------------------------------------------------------------
c
      data tcfil /'TOPCONC   '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- set the string for error messages ---
c
      action = 'Reading the NetCDF top boundary coditions file.'
c
c --- get the type of file to make sure it is a boundary file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(itc, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(itc, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7002
      endif
      if( name_in(:7) .NE. tcfil(:7) ) goto 7001
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(itc,action,1,.FALSE.,.TRUE.,.FALSE.,.FALSE.,itmp)
c
c --- call routine to make sire file spans the episode ---
c
      call ncf_chk_tstep(itc,action,begdate,begtim,enddate,endtim,.FALSE.)
c
c --- call routine to setup species mappping array ---
c
      call ncf_set_species_mapping(itc,action,spname,nspec,ntcspc,ltcmap)
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_BNDPREP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get necessary global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_BNDPREP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',tcfil(:istrln(tcfil))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in NCF_BNDPREP:'
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
