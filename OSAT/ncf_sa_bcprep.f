      subroutine ncf_sa_bcprep(begtim,begdate,endtim,enddate)
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
c     NCF_SA_BCPREP reads and checks the SA boundary condition  file.
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
      character*(16*MXTRSP) varlist
      character*200         action, fname
      character*10          bcfiltyp, name_in, this_var
      character*3           first_class
      integer               this_dimid, ierr, this_varid, itmp, nvars_in, l
c
      character*10,     allocatable :: bcspec(:)
c
c-----------------------------------------------------------------------
c    Data statments:
c-----------------------------------------------------------------------
c
      data bcfiltyp /'BOUNDARY  '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- skip if this is a NetCDF file ---
c
      action = 'Reading NetCDF SA boundary conditions file.'
c
c --- get the type of file to make sure it is a point source file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iorbc, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iorbc, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7001
      endif
      if( name_in(:8) .NE. bcfiltyp(:8) ) goto 7000
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(iorbc,action,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,itmp)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iorbc,action,begdate,begtim,
     &                                        enddate,endtim,.FALSE.)
c
c --- get the number of variables in the file ----
c
      this_var = 'NVARS'
      ierr = nf_get_att_int(iorbc, NF_GLOBAL, this_var, nvars_in)
      if( ierr .NE. NF_NOERR ) goto 7002
c
c  --- check that this number of tracer species is allowed ---
c
      if( nvars_in .GT. MXTRSP ) then
         write(iout,'(//,A)') 'ERROR in NCF_SA_BCPREP:'
         write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
         write(iout,*) 'Please change the value for parameter: MXTRSP'
         write(iout,*) 'It should be set to a value of at least: ',nvars_in
         call flush(iout)
         call camxerr()
      endif
c
c  --- allocate and read the list of variable names ----
c
      varlist = ' '
      this_var = 'VAR-LIST'
      ierr = nf_get_att_text(iorbc, NF_GLOBAL, this_var, varlist)
      if( ierr .NE. NF_NOERR ) goto 7002
      num_iorbc = 0
      do l=1,nvars_in
         if( varlist(1+(l-1)*16:16+(l-1)*16) .NE. '                ' ) 
     &                                         num_iorbc = num_iorbc + 1
      enddo
      allocate( spc_iorbc(num_iorbc) )
      first_class = varlist(1:3)
      do l=1,num_iorbc
         spc_iorbc(l) = varlist(1+(l-1)*16:10+(l-1)*16)
         if( spc_iorbc(l)(1:3) .EQ. first_class ) 
     &                                         ncls_iorbc = ncls_iorbc + 1
      enddo
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_SA_BCPREP: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',bcfiltyp(:istrln(bcfiltyp))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_SA_BCPREP: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)',ERR=9999) 'Cannot find global variable for ',
     &                                       'type of file: CAMx_NAME'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_SA_BCPREP: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)',ERR=9999) 'Cannot find global variable: ',
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
