      subroutine nodes_tstep(numprocs,iproc_id)
c
c----CAMx v7.20 220430
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

      implicit none
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
c        iproc_id            process number for this process
c     Output:
c
c    Called by:
c       CAMX
c    Subroutines called:
c       NODE_RECV_1SPECIES_DATA
c       MASTER_SEND_1SPECIES_DATA
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c       03/15/09     Added code for deposition output for tracers
c       10/29/09     Added code for RTRAC surface model
c       04/30/13     Added surface model
c       08/26/21     RTRAC surface file doubles as (default)
c                    dep output or (flagged) surface model output
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'camx_aero.inc'
      include 'soap.inc'
      include 'rtracsrf.inc'
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
      integer :: ilay
      integer :: ispc
      integer :: i
      integer :: nlayav
      integer :: iptr_rtsolmas
c    
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      iptr_rtsolmas = 1
      if( .NOT. lmpi ) return
c
c  --- pass the fields that are dimensioned by mvecem ---
c
      do i=1,ngrid      
         do ilay=1,nlayers_ems
            do ispc=1,nemspc
               if( iproc_id .gt. 0 ) then
                 call node_recv_1species_data  (aremis(iptrem(i)), i,
     &                             nlayers_ems, ilay, nemspc, ispc, itag)
               else
                 call master_send_1species_data(aremis(iptrem(i)), i,
     &                             nlayers_ems, ilay, nemspc, ispc, itag)
               endif
               itag = itag + 1
            enddo
         enddo
c
c  --- pass the concentration fields ----
c
         do ilay=1,nlay(i)
            do ispc=1,nspec
               if( iproc_id .gt. 0 ) then
                 call node_recv_1species_data  (conc(iptr4d(i)), i,
     &                              nlay(i), ilay, nspec, ispc, itag)
               else
                 call master_send_1species_data(conc(iptr4d(i)), i,
     &                             nlay(i), ilay, nspec, ispc, itag)
               endif
               itag = itag + 1
            enddo
         enddo
c
         nlayav = nlay(i)
         if( .NOT. l3davg(i) ) nlayav = 1
         do ilay=1,nlayav
            do ispc=1,navspc
                if( iproc_id .gt. 0 ) then
                  call node_recv_1species_data  (avcnc(iptrav(i)), i,
     &                               nlayav, ilay, navspc, ispc, itag)
                else
                  call master_send_1species_data(avcnc(iptrav(i)), i,
     &                               nlayav, ilay, navspc, ispc, itag)
                endif
                itag = itag + 1
            enddo
         enddo
c
c
c  --- pass the surface model soil mass fields ----
c
         if (lsrfmod) then
           do ispc=1,nsmspc
              if( iproc_id .gt. 0 ) then
                call node_recv_1species_data (solmas(iptrsm(i)), i,
     &                                       1, 1, nsmspc, ispc, itag)
              else
                call master_send_1species_data(solmas(iptrsm(i)), i,
     &                                       1, 1, nsmspc, ispc, itag)
              endif
              itag = itag + 1
           enddo
c
c  --- pass the surface model veg mass fields ----
c
           do ispc=1,nsmspc
              if( iproc_id .gt. 0 ) then
                call node_recv_1species_data (vegmas(iptrsm(i)), i,
     &                                       1, 1, nsmspc, ispc, itag)
              else
                call master_send_1species_data(vegmas(iptrsm(i)), i,
     &                                       1, 1, nsmspc, ispc, itag)
              endif
              itag = itag + 1
           enddo
         endif
c
c  --- pass the emissions array for probing tools ----
c
         if( ltrace .OR. lddm .OR. lhddm ) then 
c
            do ilay=1,nlayers_ems
               do ispc=1,ntotsp
                  if( iproc_id .gt. 0 ) then
                    call node_recv_1species_data(saemis(ipsa3d_ems(i)), i,
     &                               nlayers_ems, ilay, ntotsp, ispc, itag)
                  else
                    call master_send_1species_data(saemis(ipsa3d_ems(i)), i,
     &                               nlayers_ems, ilay, ntotsp, ispc, itag)
                  endif
                  itag = itag + 1
               enddo
            enddo
c
c  --- pass the concentration fields for probing tools ----
c
            do ilay=1,nlay(i)
               do ispc=1,ntotsp
                   if( iproc_id .gt. 0 ) then
                     call node_recv_1species_data  (ptconc(ipsa3d(i)), i,
     &                                  nlay(i), ilay, ntotsp, ispc, itag)
                   else
                     call master_send_1species_data(ptconc(ipsa3d(i)), i,
     &                                 nlay(i), ilay, ntotsp, ispc, itag)
                   endif
                   itag = itag + 1
               enddo
            enddo
c
c  --- pass the average concentration fields for probing tools ----
c
            nlayav = nlay(i)
            if( .NOT. lsa_3davrg ) nlayav = 1
            do ilay=1,nlayav
               do ispc=1,ntotsp
                   if( iproc_id .gt. 0 ) then
                        call node_recv_1species_data  (ptavrg(ipsa2d_avrg(i)), i,
     &                                     nlayav, ilay, ntotsp, ispc, itag)
                   else
                     call master_send_1species_data(ptavrg(ipsa2d_avrg(i)), i,
     &                                     nlayav, ilay, ntotsp, ispc, itag)
                   endif
                   itag = itag + 1
               enddo
            enddo
c
c  --- pass the dry depostion fields for probing tools ----
c
            if( lptdepout ) then
              do ispc=1,notimespc
                  if( iproc_id .gt. 0 ) then
                     call node_recv_1species_data  (ptdryfld(ipsadep(i)), i,
     &                                           1, 1, notimespc, ispc, itag)
                  else
                    call master_send_1species_data(ptdryfld(ipsadep(i)), i,
     &                                           1, 1, notimespc, ispc, itag)
                  endif
                  itag = itag + 1
              enddo
c
c  --- pass the wet depostion fields for probing tools ----
c
              do ispc=1,notimespc
                  if( iproc_id .gt. 0 ) then
                     call node_recv_1species_data  (ptwetfld(ipsadep(i)), i,
     &                                           1, 1, notimespc, ispc, itag)
                  else
                    call master_send_1species_data(ptwetfld(ipsadep(i)), i,
     &                                           1, 1, notimespc, ispc, itag)
                  endif
                  itag = itag + 1
              enddo
            endif
c
c  --- pass the RTRAC soil mass fields ----
c
            if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
              do ispc=1,ntotsp
                  if( iproc_id .gt. 0 ) then
                     call node_recv_1species_data (rtsolmas(iptr_rtsolmas), i,
     &                                           1, 1, ntotsp, ispc, itag)
                  else
                    call master_send_1species_data(rtsolmas(iptr_rtsolmas), i,
     &                                           1, 1, ntotsp, ispc, itag)
                  endif
                  itag = itag + 1
              enddo
c
c  --- pass the RTRAC veg mass fields ----
c
              do ispc=1,ntotsp
                  if( iproc_id .gt. 0 ) then
                     call node_recv_1species_data (rtvegmas(iptr_rtsolmas), i,
     &                                           1, 1, ntotsp, ispc, itag)
                  else
                    call master_send_1species_data(rtvegmas(iptr_rtsolmas), i,
     &                                           1, 1, ntotsp, ispc, itag)
                  endif
                  itag = itag + 1
              enddo
            endif
         endif
c
c  --- pass the IRR fields for process analysis ---
c
         if( lirr ) then
            do ilay=1,nlay(i)
              do ispc=1,ntotsp
                   if( iproc_id .gt. 0 ) then
                     call node_recv_1species_data  (ptconc(ipsa3d(i)), i,
     &                                  nlay(i), ilay, ntotsp, ispc, itag)
                   else
                     call master_send_1species_data(ptconc(ipsa3d(i)), i,
     &                                 nlay(i), ilay, ntotsp, ispc, itag)
                   endif
                   itag = itag + 1
               enddo
            enddo
         endif
      enddo
c
c  --- return to the calling routine ---
c
      end
