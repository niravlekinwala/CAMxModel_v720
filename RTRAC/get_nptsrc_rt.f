      subroutine get_nptsrc_rt( )
      use filunit
      use ptemiss
      use tracer
      implicit none
c 
c----CAMx v7.20 220430
c 
c     GET_NPTSRC_RT loops over all of the elevated point source files 
c     for the source groups and calculates the total number of additonal
c     sources. This is to allocate the point source arrays appropriately.
c                           
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c 
c     Input arguments: 
c             
c     Output arguments: 
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        PNTPREP
c 
      include 'camx.prm'
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----External functions
c
      integer istrln
c
c-----Local variables
c
      character*300 fname, action
      character*10  this_var
      integer       numpts, igroup, idxfile, num_ptsfiles, iounit
      integer       irec, ipt, ierr, this_dimid, this_varid
      integer       idxpt, idum
      real          tdummy
      logical       lbinary
c
c-----Entry point
c
c-----Allocate array for list of all points
c
c-----Loop over groups, group 0 is regular model
c
      numpts = 0
c
c-----Loop over groups, group 0 is regular model
c
      do igroup=1,ngroup
c
c-----Loop over number of files for this group
c
         num_ptsfiles = num_iortpt(igroup)
         do idxfile=1,num_ptsfiles
c
c-----Only continue if point source file is supplied
c
            if( .NOT. lptsrc ) goto 9999
c
c-----figure out how to deal with this file
c
            lbinary = .TRUE.
            if( .NOT. ltptfl(igroup,idxfile) ) cycle
            if( is_netcdf_iortpt(igroup,idxfile) ) lbinary = .FALSE.
            iounit = iortpt(igroup,idxfile)
            fname = tptfil(igroup,idxfile)
            write(action,'(2A)') 'Cannot read header of point source file: ',
     &                                                  fname(:istrln(fname))
c
c----process the binary file
c
            if( lbinary ) then
               rewind(iounit)
               do irec=1,4
                  read(iounit,ERR=7000)
               enddo
               read(iounit,ERR=7000) idum,numpts
               nptsrc_safile(igroup,idxfile) = numpts
               nptsrc = nptsrc + numpts
               if( nptsrc .GT.  MXPTSRC ) goto 7002
c
c-----this file is NetCDF format ---
c
            else
              ierr = nf_inq_dimid(iounit, "COL", this_dimid )
              if( ierr .NE. NF_NOERR ) goto 7001
              ierr = nf_inq_dimlen(iounit,this_dimid,numpts)
              if( ierr .NE. NF_NOERR ) goto 7001
              nptsrc_safile(igroup,idxfile) = numpts
              nptsrc = nptsrc + numpts
              if( nptsrc .GT.  MXPTSRC ) goto 7002
            endif
c
c-----if regular model inventory file just add to list ----
c
            nptsrc_safile(igroup,idxfile) = numpts
            if( numpts .GT.  MXPTSRC ) goto 7002
c
c----next file in this group
c
         enddo
c
c----next group
c
      enddo
c
      goto 9999
c
c---- error messages
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in GET_NPTSRC_RT:'
      write(iout,'(1X,2A)') action
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in GET_NPTSRC_RT'
      write(iout,'(A)') action
      write(iout,'(A)') 'Cannot find the dimension id for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in GET_NPTSRC_RT:'
      write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
      write(iout,*) 'Please change the value for parameter: ','MXPTSRC'
      write(iout,*) 'It should be set to a value of at least: ',nptsrc
      call flush(iout)
      call camxerr()
c
 9999 continue
      return
      end
