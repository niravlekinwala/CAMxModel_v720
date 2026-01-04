      subroutine ncf_rdpthdr_sa(begtim,begdate,endtim,enddate)
      use filunit
      use ptemiss
      use tracer
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_RDPTHDR_SA reads the header information from the global section
c     of the NetCDF point source files for SA and checks for consistency 
c     with model simulation.
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
      character*10  pntfil, name_in, this_var
      integer       this_dimid, idxfile, ierr, this_varid, tmpmap(MXSPEC)
      integer       igroup, iounit, numpts
c
c-----------------------------------------------------------------------
c    Data statments:
c-----------------------------------------------------------------------
c
      data pntfil /'PTSOURCE  '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      numpts = 0
c
c   --- loop over groups and number of files ---
c
      do igroup=1,ngroup
      do idxfile=1,num_iortpt(igroup)
        iounit = iortpt(igroup,idxfile)
c
c   --- skip if this is a NetCDF file ---
c
         if( .NOT. is_netcdf_iortpt(igroup,idxfile) ) cycle
         action = 'Reading NetCDF point source file.'
c
c --- get the type of file to make sure it is a point source file ---
c
         this_var = 'CAMx_NAME'
         ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) then
            this_var = 'NAME'
            ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
            if( ierr .NE. NF_NOERR ) goto 7007
         endif
         if( name_in(:8) .NE. pntfil(:8) ) goto 7001
c
c --- call routine to make sure file spans the episode ---
c
         call ncf_chk_tstep(iounit,action,begdate,begtim,
     &                                        enddate,endtim,le1day)
c
c --- call routine to setup species mappping array ---
c
         call ncf_set_emiss_mapping(iounit,action)
c
c --- get the number of points (used as dimension for columns) ---
c
         ierr = nf_inq_dimid(iounit, "COL", this_dimid )
         if( ierr .NE. NF_NOERR ) goto 7002
         ierr = nf_inq_dimlen(iounit,this_dimid, 
     &                              nptsrc_safile(igroup,idxfile))
         if( ierr .NE. NF_NOERR ) goto 7002
c
c --- if this is the first time, set the conters ---
c
         idx_start_sapts(igroup,idxfile) = numpts
         numpts = numpts + nptsrc_safile(igroup,idxfile)
c
      enddo
      enddo
c
c-----------------------------------------------------------------------
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get necessary global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',pntfil(:istrln(pntfil))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension value for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7007 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDPTHDR:'
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
