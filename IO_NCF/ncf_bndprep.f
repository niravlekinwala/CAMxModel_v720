c**** NCF_BNDPREP
c
      subroutine ncf_bndprep(begtim,begdate,endtim,enddate)
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
c     NCF_BNDPREP reads the boundary conditions NETCDF file and checks
c     the global varibles for consistency with user defined values.
c     It also prepares mapping variables that map the BC species list
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
      character*10  bcfil, name_in, this_var
      integer       ierr, i, j, n, itmp
c
c-----------------------------------------------------------------------
c    Data Statements:
c-----------------------------------------------------------------------
c
      data bcfil /'BOUNDARY  '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- set the string for error messages ---
c
      action = 'Reading the NetCDF lateral boundary coditions file.'
c
c --- get the type of file to make sure it is a boundary file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(ibc, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(ibc, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7002
      endif
      if( name_in(:8) .NE. bcfil(:8) ) goto 7001
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(ibc,action,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,itmp)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(ibc,action,begdate,begtim,enddate,endtim,.FALSE.)
c
c --- call routine to setup species mappping array ---
c
      call ncf_set_species_mapping(ibc,action,spname,nspec,nbcspc,lbcmap)
c
c --- set the values for boundary cell definition to edge cells ---
c    
      do j=1,nrow(1)
         ibeg(j) = 2
         iend(j) = ncol(1)-1
      enddo
      do i=1,ncol(1)
         jbeg(i) = 2
         jend(i) = nrow(1)-1
      enddo
c
c --- Ensure that rows 1 & ny, and columns 1 & nx are boundary cells ---
c
      ibeg(1) = -999
      iend(1) = -999
      ibeg(nrow(1)) = -999
      iend(nrow(1)) = -999
      jbeg(1) = -999
      jend(1) = -999
      jbeg(ncol(1)) = -999
      jend(ncol(1)) = -999
c
c --- echo coarse grid boundary definition, and check that if nested grids
c     are specified, their boundaries are not outside the coarse
c     grid region
c
      do j = 1,nrow(1)
        do n = 2,ngrid
          if (j.lt.jnst1(n) .or. j.gt.jnst2(n)) cycle
          if (ibeg(j).eq.-999 .or. inst1(n).lt.ibeg(j) .or.
     &        inst2(n).gt.iend(j)) then
            write(iout,'(//,a)') 'ERROR in NCF_BNDPREP:'
            write(iout,*)'Nested grid extends into coarse grid',
     &                   ' boundary region -- Nest = ',n
            write(iout,*)'j,ibeg,iend,inst1,inst2: ',
     &                    j,ibeg(j),iend(j),inst1(n),inst2(n)
            call camxerr()
          endif
        enddo
      enddo
      do i = 1,ncol(1)
        do n = 2,ngrid
          if (i.lt.inst1(n) .or. i.gt.inst2(n)) cycle
          if (jbeg(i).eq.-999 .or. jnst1(n).lt.jbeg(i) .or.
     &        jnst2(n).gt.jend(i)) then
            write(iout,'(//,a)') 'ERROR in NCF_BNDPREP:'
            write(iout,*)'Nested grid extends into coarse grid',
     &                   ' boundary region -- Nest = ',n
            write(iout,*)'i,jbeg,jend,jnst1,jnst2: ',
     &                    i,jbeg(i),jend(i),jnst1(n),jnst2(n)
            call camxerr()
          endif
        enddo
      enddo
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
      write(iout,'(2A)') 'Looking for type: ',bcfil(:istrln(bcfil))
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
