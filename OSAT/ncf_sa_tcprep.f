      subroutine ncf_sa_tcprep(begtim,begdate,endtim,enddate)
      use filunit
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_SA_TCPREP reads and checks the SA boundary condition  file.
c     It will determine the species in the file and populate the 
c     tracer array appropriately.
c
c      Copyright 1996 - 2022
c     Ramboll
c
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
      include 'flags.inc'
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
      character*200 action, fname
      character*10  tcfiltyp, name_in, this_var
      character*3   first_class
      integer       this_dimid, ierr, this_varid, itmp, l
c
      character*10, allocatable :: tcspec(:)
c
c-----------------------------------------------------------------------
c    Data statments:
c-----------------------------------------------------------------------
c
      data tcfiltyp /'TOPCONC   '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- skip if this is a NetCDF file ---
c
      action = 'Reading NetCDF SA top boundary conditions file.'
c
c --- get the type of file to make sure it is a point source file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iortc, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iortc, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7001
      endif
      if( name_in(:8) .NE. tcfiltyp(:8) ) goto 7000
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(iortc,action,1,.FALSE.,.TRUE.,.FALSE.,.FALSE.,itmp)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iortc,action,begdate,begtim,
     &                                        enddate,endtim,.FALSE.)
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_SA_TCPREP: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',tcfiltyp(:istrln(tcfiltyp))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_SA_TCPREP: '
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
