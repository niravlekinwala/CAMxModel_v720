      subroutine get_nptsrc( )
      use filunit
      use chmstry
      use ptemiss
      implicit none
c 
c----CAMx v7.20 220430
c 
c     GET_NPTSRC loops over all of the elevated point source files 
c     and retrieves the number of point sources in each file. It then 
c     calculates the total number in the inventory.
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
      include 'netcdf.inc'
c
c-----Local variables
c
      integer idxfile, ierr, this_dimid, irec, idum
c
c-----Entry point
c
c-----Loop over all of the files
c
      do idxfile=1,npoint_files
c
c-----this file is old style binary ---
c
        if( .NOT. is_netcdf_iptem(idxfile) ) then
c
c-----Read 1st PT header record and check inputs 
c             
           rewind(iptem(idxfile))
           do irec=1,4
              read(iptem(idxfile))
           enddo
c 
           read(iptem(idxfile)) idum,nptsrc_files(idxfile)
        else
c
c-----this file is NetCDF format ---
c
           ierr = nf_inq_dimid(iptem(idxfile), "COL", this_dimid )
           if( ierr .NE. NF_NOERR ) goto 7000
           ierr = nf_inq_dimlen(iptem(idxfile),this_dimid,
     &                                    nptsrc_files(idxfile))
           if( ierr .NE. NF_NOERR ) goto 7000
        endif
      enddo
c
c---- all files processed, now calculate starting points ---
c
      idx_start_pts(1) = 0
      nptsrc = nptsrc_files(1)
      do idxfile=2,npoint_files
         idx_start_pts(idxfile) = nptsrc
         nptsrc = nptsrc + nptsrc_files(idxfile)
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
