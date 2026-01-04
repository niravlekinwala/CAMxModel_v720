c***** NCF_EMPREPDDM
c
      subroutine ncf_emprepddm(igrid,begdate,begtim,enddate,endtim)
      use filunit
      use tracer
      use ptemiss
      use chmstry
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
c   This version is for DDM.
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
      integer       num_emsfiles, num_ptsfiles, igroup, idxfile
      integer       numpts_tot
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      numpts_tot = idx_start_pts(npoint_files) + nptsrc_files(npoint_files)
c
c   --- initialize postion zero, which is for the model emissions ----
c
      do igroup=1,ngroup
          if( larsrc ) then
             num_emsfiles = num_iortem(igrid,igroup)
             if( igroup .EQ. 0 ) num_emsfiles = nemiss_files(igrid)
             do idxfile=1,num_emsfiles
                if( is_netcdf_iortem(igrid,igroup,idxfile)) then
                   write(action,'(2A,I3,A,I3,A,I3)') 'Reading the NetCDF ',
     &                        'gridded DDM emissions file. Grid: ',
     &                        igrid,' Gourp: ',igroup,' File: ',idxfile
                   call ncf_areaprep(action,iortem(igrid,igroup,idxfile),
     &                                igrid,begtim,begdate,endtim,enddate,
     &                                 buffer_offset_iortem(igrid,igroup,idxfile))
                endif
             enddo
          endif
c
c   --- point sources are only for the coarse grid -- grid #1 ---
c
          if( lptsrc  .AND. igrid .EQ. 1 ) then
             num_ptsfiles = num_iortpt(igroup)
             do idxfile=1,num_ptsfiles
                if( is_netcdf_iortpt(igroup,idxfile) )
     &              call ncf_rdpthdr_ddm(igroup,idxfile,begtim,begdate,
     &                                        endtim,enddate,numpts_tot)
             enddo
          endif
c
      enddo
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
