      subroutine nodes_topc(numprocs,iproc_id)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use grid
      use chmstry
      use filunit
      use o3colmap
      use bndary
      use camxfld
      use camxcom
      use pigsty
      use ptemiss
      use procan
      use rtracchm
      use tracer
c
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c        This routine passes all of the data that is time-step
c        dependent to the compute nodes when in MPI mode.
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c        numprocs            the number of processes
c        iproc_id            processes number of current process
c     Output:  
c
c    Called by:
c       BNDRY_UPDT
c    Subroutines called:
c       NODE_RECV_1SPECIES_DATA
c       MASTER_SEND_1SPECIES_DATA
c       NODES_PASS
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     11/28/16   --bkoo--       Updated DDM for new top con
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'camx_aero.inc'
      include 'soap.inc'
      include 'mpif.h'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer :: numprocs
      integer :: iproc_id
c     
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: mvec2d
      integer :: mvecem
      integer :: ispc
      integer :: i
c     
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if (.not. lmpi .and. iproc_id .gt. 0) return
c
c   --- calculate the number of bytes to send ---
c
      mvec2d = 0
      do i=1,ngrid
         mvec2d = mvec2d + nrow(i) * ncol(i)
      enddo
      mvecem = mvec2d * nspec
c
c  --- pass the fields that are dimensioned by mvecem ---
c
      do i=1,ngrid
         do ispc=1,nspec
            if (iproc_id .gt. 0) then
               call node_recv_1species_data  (ctop(iptr1lay(i)),i,
     &                                        1,1,nspec,ispc,itag )
            else
               call master_send_1species_data(ctop(iptr1lay(i)),i,
     &                                        1,1,nspec,ispc,itag )
            endif
            itag = itag + 1
         enddo
      enddo
c
c  --- If doing the probing tools, pass those fields ---
c
      if ( lddm .OR. lhddm ) then
         do i=1,ngrid
            do ispc=1,ntotsp
               if (iproc_id .gt. 0) then
                  call node_recv_1species_data  (pttop(ipsa2d(i)),i,
     &                                           1,1,ntotsp,ispc,itag)
               else
                  call master_send_1species_data(pttop(ipsa2d(i)),i,
     &                                           1,1,ntotsp,ispc,itag)
               endif
               itag = itag + 1
            enddo
         enddo
      endif
c
      end
