c***** NCF_EMPREPSA
c
      subroutine ncf_emprepsa(igrid,begdate,begtim,enddate,endtim)
      use filunit
      use chmstry
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the NetCDF source group emissions files and
c   verifies that all files are consistent with the model simulation.
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Inputs:
c         igrid   I grid number
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
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      real    begtim
      integer begdate
      real    endtim
      integer enddate
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      integer       num_emsfiles, num_ptsfiles, igroup, idxfile, i
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- initialize postion zero, which is for the model emissions ----
c
      if( larsrc ) then
         do igroup=1,ngroup
             num_emsfiles = num_iortem(igrid,igroup)
             do idxfile=1,num_emsfiles
                if( is_netcdf_iortem(igrid,igroup,idxfile)) then
                   write(action,'(2A,I3,A,I3,A,I3)') 'Reading the NetCDF ',
     &                        'gridded SA emissions file. Grid: ',
     &                        igrid,' Group: ',igroup,' File: ',idxfile
                   call ncf_areaprep(action,iortem(igrid,igroup,idxfile),
     &                               igrid,begtim,begdate,endtim,enddate,
     &                                   buffer_offset_iortem(igrid,igroup,idxfile))
                endif
             enddo
         enddo
      endif
c
c  --- update index for group zero ---
c
      if( lptsrc ) then
         do idxfile=1,npoint_files
            do i=1,nptsrc_files(idxfile)
               idx_point_in_list(0,idxfile,i) = i + idx_start_pts(idxfile)
            enddo
         enddo
      endif
c
c  ---- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
