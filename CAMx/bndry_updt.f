      subroutine bndry_updt(numcols,numrows,numlays,bndtim,bnddate,
     &                      toptim,topdate,nsteps,numprocs,iproc_id)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use filunit
      use camxcom
      use camxfld
      use chmstry
      use grid
      use bndary
      use node_mod
      use tracer
c
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c
c    Called by:
c       CAMX
c    Subroutines called:
c       READBND
c       READTOP
c       FLUSH
c       NODE_RECV_1SPECIES_DATA
c       MASTER_SEND_1SPECIES_DATA
c       RDBCRT
c       NCF_RDBCRT
c       CLRBDYSA
c       FILBDYSA
c       CLRBDYDDM
c       ZEROS
c       RDBCDDM
c       NODES_PASS
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     11/4/09 -cemery- Removed input top concentrations
c      4/7/14 -cemery- Added top con file
c     11/28/16 -bkoo-  Updated DDM for new top con
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'mpif.h'
c
c========================= Probing Tool End ============================
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer :: numcols(*)
      integer :: numrows(*)
      integer :: numlays(*)
c
      real    :: bndtim
      real    :: toptim
c
      integer :: bnddate
      integer :: topdate
      integer :: nsteps
      integer :: numprocs
      integer :: iproc_id
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer (kind=8) :: mvsa3d
      integer          :: ncola
      integer          :: nrowa
      integer          :: nmx1d
      integer          :: nedge
      integer          :: i
      integer (kind=8) :: idxsa
      integer          :: ierr
      integer (kind=8) :: numsa
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- If reading NetCDF BC file call routine to see if update is needed ---
c
      if( is_netcdf_ibc .AND. iproc_id .EQ. 0) 
     &                   call ncf_readbnd(time,date,did_update_ibc)
      if( lmpi .AND. is_netcdf_ibc ) then
        call nodes_pass(did_update_ibc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
        if( did_update_ibc ) then
           nedge = MAX(ncol(1),nrow(1))
           if( iproc_id .EQ. 0 ) 
     &                call master_bndconc(ncol(1),nrow(1),nlay(1),
     &                                          nspec,nedge,conc,bndconc)
            call edge_pass(nodeedge,bndconc,4*nedge*nlay(1)*nspec,MPI_REAL,
     &                                             itag,numprocs,iproc_id) 
            if( iproc_id .NE. 0 ) then
                if( nodeedge(iproc_id) ) 
     &                     call nodes_bndconc(mi0(1),mj0(1),
     &                        mmxp(1),mmyp(1),ncol(1),nrow(1),nlay(1),
     &                                         nspec,nedge,conc,bndconc)
           endif
        endif
      endif
c
c  --- Check if master grid boundary data are to be read ---
c
      if (.NOT. is_netcdf_ibc .AND. date .EQ. bnddate 
     &                           .AND. ABS(time-bndtim) .LT. 0.01) then
         if (iproc_id .EQ. 0) then
c
c========================= Source Apportion End ======================
c
            write(*,'(a20,$)') 'readbnd ......'
            call readbnd(bndtim,bnddate)
            write(*,'(a)') '   Done'
            call flush(6)
         endif
c
c   --- if doing MPI, update edge cells on each slice ---
c
         if( lmpi ) then
           nedge = MAX(ncol(1),nrow(1))
           if( iproc_id .EQ. 0 ) 
     &                call master_bndconc(ncol(1),nrow(1),nlay(1),
     &                                          nspec,nedge,conc,bndconc)
            call edge_pass(nodeedge,bndconc,4*nedge*nlay(1)*nspec,MPI_REAL,
     &                                             itag,numprocs,iproc_id) 
            if( iproc_id .NE. 0 ) then
                   if( nodeedge(iproc_id) ) 
     &                     call nodes_bndconc(mi0(1),mj0(1),
     &                        mmxp(1),mmyp(1),ncol(1),nrow(1),nlay(1),
     &                                         nspec,nedge,conc,bndconc)
           endif
         endif
c
c========================= Source Apportion Begin ======================
c
c  --- call routine to clear the old boundary cells and then
c      call routine to fill with new boundary concentrations ---
c
c  --- For the Master Grid ---
c
         if (iproc_id .eq. 0) then
            if (ltrace) then
               if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
                  call rdbcrt(nsteps,numcols(1),numrows(1),numlays(1),
     &                        ntotsp,ptconc(1),tempk(1),press(1)      )
                  call ncf_rdbcrt(numcols(1),numrows(1),numlays(1),
     &                        ntotsp,ptconc(1),tempk(1),press(1)      )
               else
                  if( .NOT. lsa_iorbc ) then
                      call clrbdysa(1,numcols(1),numrows(1),
     &                          numlays(1),ntotsp,ptconc(1))
                      call filbdysa(1,numcols(1),numrows(1),numlays(1),
     &                          nspec,ntotsp,conc(1),ptconc(1)     )
                  endif
               endif
c
c  --- make sure it is not below lower bound ---
c
               mvsa3d = 0
               ncola = 0
               nrowa = 0
               do i=1,ngrid
                  ncola = MAX(ncola,numcols(i))
                  nrowa = MAX(nrowa,numrows(i))
                  mvsa3d = mvsa3d + DBLE(numrows(i)) * DBLE(numcols(i)) * DBLE(numlays(i))
               enddo
               nmx1d  = MAX( nrowa, ncola )
               numsa = mvsa3d*DBLE(ntotsp)
               do idxsa=1,numsa
                  ptconc(idxsa) = AMAX1(ptconc(idxsa),BNDLPT)
               enddo
            endif
c
c========================= Source Apportion End ========================
c
c
c============================= DDM Begin ===============================
c
c  --- call routine to clear the old boundary cells and then
c      call routine to fill with new boundary concentrations ---
c
c  --- For the Master Grid ---
c
            if( (lddm.OR.lhddm) .AND. nbcddm .GT. 0 ) then
               call rdbcddm(numcols(1),numrows(1),numlays(1),
     &                      ntotsp,ptconc(1),tempk(1),press(1))
               call ncf_rdbcddm(numcols(1),numrows(1),numlays(1),
     &                      ntotsp,ptconc(1),tempk(1),press(1))
            endif
         endif
c
c  --- update the DDM boundary sensitivites for the MPI slices ---
c
         if( (lddm.OR.lhddm) .AND. lmpi ) then
             nedge = MAX(ncol(1),nrow(1))
             if( iproc_id .EQ. 0 )
     &                call master_bndconc(ncol(1),nrow(1),nlay(1),
     &                                    ntotsp,nedge,ptconc,bndddm)
             call edge_pass(nodeedge,bndddm,4*nedge*nlay(1)*ntotsp,
     &                               MPI_REAL,itag,numprocs,iproc_id)
             if( iproc_id .NE. 0 ) then
               if( nodeedge(iproc_id) )
     &                call nodes_bndconc(mi0(1),mj0(1),
     &                      mmxp(1),mmyp(1),ncol(1),nrow(1),nlay(1),
     &                                     ntotsp,nedge,ptconc,bndddm)
             endif
         endif
c
c============================= DDM End =================================
c
c========================= Source Apportion Begin ======================
c
c  --- update the RTRAC boundary conditions for the MPI slices ---
c
         if( lmpi .AND. ltrace .AND.
     &           (tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) ) then
             nedge = MAX(ncol(1),nrow(1))
             if( iproc_id .EQ. 0 )
     &                call master_bndconc(ncol(1),nrow(1),nlay(1),
     &                                    ntotsp,nedge,ptconc,bndrt)
             call edge_pass(nodeedge,bndrt,4*nedge*nlay(1)*ntotsp,
     &                               MPI_REAL,itag,numprocs,iproc_id)
             if( iproc_id .NE. 0 ) then
               if( nodeedge(iproc_id) )
     &                call nodes_bndconc(mi0(1),mj0(1),
     &                      mmxp(1),mmyp(1),ncol(1),nrow(1),nlay(1),
     &                                     ntotsp,nedge,ptconc,bndrt)
             endif
         endif
c
c========================= Source Apportion End ========================
c
      endif
c
c========================= Source Apportion Begin ======================
c
      if( ltrace .AND. lsa_iorbc ) then
         if( iproc_id .EQ. 0 ) then
            if( is_netcdf_iorbc ) then
               call ncf_rd_sa_bcfile(bndtim,bnddate)
            else
               if( date .EQ. sa_bnddate
     &                      .AND. ABS(time-sa_bndtim) .LT. 0.01) then
                  call rd_sa_bcfile(bndtim,bnddate)
               endif
            endif
            if( lsa_iortc ) then
               if( lsa_iortc ) then
                  call ncf_rd_sa_tcfile(bndtim,bnddate)
               else
                   if( date .EQ. sa_bnddate
     &                      .AND. ABS(time-sa_bndtim) .LT. 0.01) then
                      call rd_sa_tcfile(bndtim,bnddate)
                   endif
               endif
            endif
         endif
         if( lmpi ) then
             nedge = MAX(ncol(1),nrow(1))
             if( iproc_id .EQ. 0 )
     &                call master_bndconc(ncol(1),nrow(1),nlay(1),
     &                                    ntotsp,nedge,ptconc,bndsa)
             call edge_pass(nodeedge,bndsa,4*nedge*nlay(1)*ntotsp,
     &                               MPI_REAL,itag,numprocs,iproc_id)
             if( iproc_id .NE. 0 ) then
               if( nodeedge(iproc_id) )
     &                call nodes_bndconc(mi0(1),mj0(1),
     &                      mmxp(1),mmyp(1),ncol(1),nrow(1),nlay(1),
     &                                     ntotsp,nedge,ptconc,bndsa)
             endif
         endif
      endif
c
c========================= Source Apportion End ======================
c
c
c  --- Check if master grid topcon data are to be read ---
c
      if (ltopcon) then
        if( is_netcdf_itc .AND. iproc_id .EQ. 0) 
     &                    call ncf_readtop(time,date,did_update_itc)
        if( lmpi .AND. is_netcdf_itc ) then
            call nodes_pass(did_update_itc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
            if( did_update_itc ) then
               call MPI_Barrier( MPI_COMM_WORLD, ierr)
               call nodes_topc(numprocs,iproc_id) ! this also updates DDM top con sens
            endif
        endif
        if( .NOT. is_netcdf_itc .AND. date .EQ. 
     &                 topdate .AND. ABS(time-toptim) .LT. 0.01) then
           if (iproc_id .EQ. 0) then
              write(*,'(a20,$)') 'readtop ......'
              call readtop(toptim,topdate)
              write(*,'(a)') '   Done'
              call flush(6)
c
c============================= DDM Begin ===============================
c
              if ((lddm.OR.lhddm) .AND. nbcddm .GT. 0) then
                call rdtcddm()
                call ncf_rdtcddm(toptim,topdate)
              endif
c
c============================= DDM End =================================
c
           endif
        endif
c
c   --- if doing MPI, update topcon on each slice ---
c
        if( lmpi ) then
          call MPI_Barrier( MPI_COMM_WORLD, ierr)
          call nodes_topc(numprocs,iproc_id) ! this also updates DDM top con sens
        endif
      endif
c
      end
