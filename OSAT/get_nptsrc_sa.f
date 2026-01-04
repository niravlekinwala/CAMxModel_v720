      subroutine get_nptsrc_sa( )
      use tracer
      use chmstry
      use ptemiss
      implicit none
c 
c----CAMx v7.20 220430
c 
c     GET_NPTSRC_SA the elevated point source files 
c     for the source groups and calculates the total number in the inventory.
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
      include 'filunit.inc'
      include 'netcdf.inc'
c
c-----Local variables
c
      integer idxfile, ierr, this_dimid, irec, idum, igroup
c
c-----Entry point
c
c-----Loop over groups
c
      do igroup=1,ngroup
c
c-----All files in this group
c
         do idxfile=1,num_iortpt(igroup)
           if( .NOT. ltptfl(igroup,idxfile) ) cycle
c
c-----this file is old style binary ---
c
           if( .NOT. is_netcdf_iortpt(igroup,idxfile) ) then
              do irec=1,4
                 read(iortpt(igroup,idxfile))
              enddo
              read(iortpt(igroup,idxfile)) idum,nptsrc_safile(igroup,idxfile)
              nptsrc = nptsrc + nptsrc_safile(igroup,idxfile)
c
c-----this file is NetCDF format ---
c
           else
              ierr = nf_inq_dimid(iortpt(igroup,idxfile), "COL", this_dimid )
              if( ierr .NE. NF_NOERR ) goto 7000
              ierr = nf_inq_dimlen(iortpt(igroup,idxfile),this_dimid,
     &                                    nptsrc_safile(igroup,idxfile))
              if( ierr .NE. NF_NOERR ) goto 7000
              nptsrc = nptsrc + nptsrc_safile(igroup,idxfile)
           endif
         enddo
      enddo
c
      goto 9999
c
c---- error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in GET_NPTSRC'
      write(iout,'(A,I3)') 'Reading NetCDF elevated point source file:',idxfile
      write(iout,'(A)') 'Cannot find the dimension id for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      return
      end
