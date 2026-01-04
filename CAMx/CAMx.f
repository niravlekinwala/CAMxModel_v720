      program CAMx
      use filunit
      use grid
      use chmstry
      use o3colmap
      use bndary
      use camxfld
      use camxcom
      use ptemiss
      use pigsty
      use procan
      use rtracchm
      use tracer
c
      use master_mod  
      use node_mod   
c
c***********************************************************************
c
c                  CCCCCCC      AA      MM      MM                       
c                 CC          AA  AA    MMM    MMM   xx   xx  
c                 CC         AA    AA   MM MMMM MM    xx xx  
c                 CC         AAAAAAAA   MM  MM  MM     xxx  
c                 CC         AA    AA   MM      MM    xx xx
c                  CCCCCCC   AA    AA   MM      MM   xx   xx
c
c        C O M P R E H E N S I V E   A I R   Q U A L I T Y   M O D E L
c        -                           -                       - 
c                         with   E X T E N S I O N S
c                                  -
c
c                          VERSION  7.20 04-30-22
c
c                            Copyright 1996 - 2022
c                                  Ramboll 
c                      7250 Redwood Blvd, Suite 105
c                             Novato, CA  94945
c                             (415) 899 - 0700
c                                www.camx.com
c
c     Revision History:
c       04/30/22   Version 7.20 released (see Release_notes_7.20)
c       01/05/21   Version 7.10 released (see Release_notes_7.10)
c       05/31/20   Version 7.00 released (see Release_notes_7.00)
c       04/30/18   Version 6.50 released (see Release_notes_6.50)
c       12/23/16   Version 6.40 released (see Release_notes_6.40)
c       04/08/16   Version 6.30 released (see Release_notes_6.30)
c       03/23/15   Version 6.20 released (see Release_notes_6.20)
c       04/02/14   Version 6.10 released (see Release_notes_6.10)
c       05/06/13   Version 6.00 released (see Release_notes_6.00)
c       11/09/12   Version 5.41 released (see Release_notes_5.41)
c       10/10/11   Version 5.40 released (see Release_notes_5.40)
c       12/01/10   Version 5.30 released (see Release_notes_5.30)
c       06/15/10   Version 5.21 released (see Release_notes_5.21)
c       06/14/10   Version 5.20.1 released (see Release_notes_5.20.1)
c       04/02/10   Version 5.20 released (see Release_notes_5.20)
c       09/09/09   Version 5.10 released (see Release_notes_5.10)
c       04/06/09   Version 5.01 released (see Release_notes_5.01)
c       05/22/08   Version 4.51 released (see Release_notes_4.51)
c       10/25/06   Version 4.40 released (see Release_notes_4.40)
c       03/29/06   Version 4.31 released (see Release_notes_4.31)
c       02/15/06   Version 4.30 released (see Release_notes_4.30)
c       06/15/05   Version 4.20 released (see Release_notes_4.20)
c       12/06/04   Version 4.11.s released (see Release_notes_4.11.s)
c       08/01/04   Version 4.10.s released (see Release_notes_4.10.s)
c       12/05/03   Version 4.03 released (see Release_notes_4.03)
c       07/09/03   Version 4.02 released (see Release_notes_4.02)
c       06/17/03   Version 4.01 released (see Release_notes_4.01)
c       05/31/03   Version 4.00 released (see Release_notes_4.00)
c       11/15/01   Version 3.10 released (see Release_notes_3.10)
c        4/25/01   Version 3.01 released (see Release_notes_3.01)
c       12/30/00   Version 3.00 released (see Release_notes_3.00) 
c       12/30/99   Version 2.03 released (see Release_notes_2.03) 
c        9/13/99   Version 2.02 released (see Release_notes_2.02) 
c        5/10/99   Version 2.01 released (see Release_notes_2.01) 
c       12/30/98   Version 2.00 released, original development
c
c***********************************************************************
c  
      include 'camx.prm'
      include 'flags.inc'
      include 'ncf_iodat.inc'
      include 'mpif.h'
      include 'rtracsrf.inc'
c
c======================== Probing Tool Begin ===========================
c
c     linit   --  logical, if true,   initialize cipr(IPR_INIT ,*,*,*)
c                          otherwise, initialize cipr(IPR_FINAL,*,*,*)
c
      logical linit
c
c========================= Probing Tool End ============================
c
      integer inpdate,emsdate,ozndate,bnddate,topdate,wrtdate,enddate 
      integer iproc_id, nlayav, ndps, num_dims
      real inptim,emstim,ozntim,bndtim,toptim,wrttim,endtim 
      character*200 action
      character*20 version
      character*8 chtime, chdate
      logical lupdtdep(MXGRID), loutspec(MXSPEC+MXTRSP)
c
      data version /'CAMx v7.20, 04-30-22'/
c
c-----Entry point
c
      write(*,'(a,a)') 'Starting ',version
      if (version .NE. PRMVERSION) then
        write(iout,'(a)') 'ERROR in CAMx:'
        write(iout,'(a)') 'CAMx version does not match CAMx.prm file.'
        write(iout,'(2a)') 'Be sure to compile with the CAMx.prm file,',
     &                   'distributed with this version of the CAMx code.'
        write(iout,'(a)') 'Do not use CAMx.prm files from older versions.'
        call camxerr()
      endif
c
c----Initialize the MPI mechanism
c
      lmpi = .TRUE.
      call MPI_INIT( ierr )
      if( ierr .NE. 0 ) then
         write(*,*) 'Note: This is not an MPI run.'
         iproc_id = 0
         numprocs = 1
      endif
      if( lmpi ) then
         call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
         call MPI_COMM_RANK( MPI_COMM_WORLD, iproc_id, ierr )
      endif
      if( numprocs .EQ. 1 ) lmpi = .FALSE.
      nsteps = 0
      itag = 1
c
c-----Start model simulation
c
      ncf_cur_tstep = 1
      if( iproc_id .EQ. 0 ) then
         write(*,*)
         call banner()
         write(*,'(a20,$)') 'Model startup ......'
      endif
      call sim_init(version,inptim,inpdate,emstim,emsdate,ozntim,
     &              ozndate,bndtim,bnddate,toptim,topdate,wrttim,
     &              wrtdate,endtim,enddate,numprocs,iproc_id) 
      nlayav = nlay(1)
      if( .NOT. l3davg(1) ) nlayav = 1
c
      if( iproc_id .EQ. 0 ) write(*,'(a)') '   Done'
      write(iout,*)
      write(idiag,*)
      write(iout,'(a)')'BEGIN SIMULATION'
      write(idiag,'(a)')'BEGIN SIMULATION'
      write(iout,'(a)') 'Starting main integration loop...'
      write(idiag,*)
      call flush(iout)
      call flush(idiag)
c
c--------------------  Main time-integration loop  ---------------------
c
  100 continue
      if( lmpi ) then
        call nodes_pass(date,1,MPI_INTEGER,itag,numprocs,iproc_id)
        call nodes_pass(time,1,MPI_REAL,itag,numprocs,iproc_id)
      endif
      call chrtime(time,date,chtime,chdate)
      if( iproc_id .EQ. 0 ) then
         nsteps = nsteps + 1
         xyordr = mod(nsteps,2)
      endif
      if( (lmpi .AND. iproc_id .EQ. 1) .OR. .NOT. lmpi ) then
         write(*,*)
         write(*,'(a,a8,1x,a8)') 'Date/time: ',chdate,chtime
         write(*,*)
      endif
c
c-----Check if time-varying grid/met data are to be read
c
      if( lmpi ) then
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        call nodes_pass(bnddate,1,MPI_INTEGER,itag,numprocs,iproc_id)
        call nodes_pass(bndtim,1,MPI_REAL,itag,numprocs,iproc_id)
        call nodes_pass(topdate,1,MPI_INTEGER,itag,numprocs,iproc_id)
        call nodes_pass(toptim,1,MPI_REAL,itag,numprocs,iproc_id)
        call nodes_pass(inpdate,1,MPI_INTEGER,itag,numprocs,iproc_id)
        call nodes_pass(inptim,1,MPI_REAL,itag,numprocs,iproc_id)
        call nodes_pass(emsdate,1,MPI_INTEGER,itag,numprocs,iproc_id)
        call nodes_pass(emstim,1,MPI_REAL,itag,numprocs,iproc_id)
        call nodes_pass(nsteps,1,MPI_REAL,itag,numprocs,iproc_id)
        call nodes_pass(xyordr,1,MPI_REAL,itag,numprocs,iproc_id)
      endif
c
      write(iout,*)
      write(iout,'(a,a8,1x,a8)') 'Date/time: ',chdate,chtime
      write(iout,*)
c
      if( date .EQ. inpdate .AND. ABS(time-inptim) .LT. 0.01) then
          do i=1,ngrid
             lupdtdep(i) = .TRUE.
          enddo
          call nodes_pass(lupdtdep,ngrid,MPI_INTEGER,
     &                                        itag,numprocs,iproc_id)
          call tstep_init(inptim,inpdate,endtim,enddate,nsteps,
     &                    numprocs,iproc_id)
      endif
c
c-----Check if emissions data are to be read
c
      if( lmpi ) call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if( date .EQ. emsdate .AND. ABS(time-emstim) .LT. 0.01)
     &                 call emiss_updt(emstim,emsdate,numprocs,iproc_id)
c
c-----Check if master grid boundary data are to be read
c
      call bndry_updt(mmxp,mmyp,mmzp,bndtim,bnddate,toptim,topdate,
     &                nsteps,numprocs,iproc_id)
c
c-----Check if ozone column data are to be read
c
      if( io3col .NE. 0 .AND. idmech .NE. 10 ) then
          call o3col_updt(ozndate,ozntim,numprocs,iproc_id)
          call nodes_pass(ozndate,1,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(ozntim,1,MPI_REAL,itag,numprocs,iproc_id)
      endif
      call newgrid(1)
c
c-----At this point we have a successful model startup
c     Stop model now if startup diagnostic flag is set
c
      if ( ldiag .and. nsteps.eq.1) then
        if( iproc_id .EQ. 0 ) then
            write(*,'(/,a,a8,1x,a8,/)')'Date/time: ',chdate,chtime
            write(*,'(a)')'SUCCESSFUL MODEL STARTUP'
            write(*,'(a)')'Diagnostic Initialization Test Complete'
        endif
        write(iout,'(/,a,a8,1x,a8,/)') 'Date/time: ',chdate,chtime
        write(iout,'(a)')'SUCCESSFUL MODEL STARTUP'
        write(iout,'(a)')'Diagnostic Initialization Test Complete'
        call flush(6)
        call MPI_Finalize(ierr)
        stop
      endif
c
c-----Update date/time for this timestep
c
      do ip=1,ngrid
        timec(ip) = time
        datec(ip) = date
      enddo
      call uptime(time,date,deltat(1))
c
c-----Initialize radical concentrations
c
      if( iproc_id .EQ. 0 .AND. lchem .and. nsteps.eq.1 
     &                                     .and. nrad.gt.0) then
        write(*,'(a20,$)') 'raddrivr ......'
        call raddrivr()
        write(*,'(a)') '   Done'
        call flush(6)
        call flush(iout)
      endif
      call chrtime(timec(1),datec(1),chtime,chdate)
c
c-----Master grid emission and transport
c
c  --- call routine to load all of the input data to
c      compute node processes ---
c
      if( lmpi .AND. nsteps .EQ. 1 ) then
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call nodes_tstep(numprocs,iproc_id)
         if( LVISPIG .AND. ipigflg .NE. 0)
     &                call nodes_met_pig(numprocs,iproc_id)
      endif
      if( .NOT.lmpi .OR. iproc_id.GT.0 ) call avgall(iproc_id,nsteps,0.)
c
      call newgrid(1)
c
      call emistrns(1,chtime,chdate,
     &                        numprocs,iproc_id,lupdtdep(1))
c
c-----Computation for nested grids
c
      if( ngrid .GT. 1 ) call nesting(numprocs,iproc_id,lupdtdep,nsteps)
c
      if( .NOT. lmpi .OR. iproc_id .GT. 0 ) then
c
         call uptime(timec(1),datec(1),deltat(1))
         call chrtime(timec(1),datec(1),chtime,chdate)
c
c-----Evolve PiG puffs on master grid
c
         if( ipigflg .NE. 0 ) call pigevol(1,iproc_id)
      endif
      if( ipigflg .NE. 0 .AND. lmpi )
     &                 call nodes_pig_pass(1,numprocs,iproc_id)
c     
c-----Master grid chemical reactions
c
      if( .NOT. lmpi .OR. iproc_id .GT. 0 ) then
         call chemrxn(1,iproc_id)
c
c-----Update average concentrations
c
         call avgall(iproc_id,nsteps,2.)
      endif
c
c-----Check if concentration fields are to be written
c
      if( lmpi ) then
         if( date .EQ. wrtdate .AND. ABS(time-wrttim) .LT. 0.01 ) then
             call MPI_Barrier(MPI_COMM_WORLD, ierr)
             call master_update(numprocs, iproc_id)
             call MPI_Barrier(MPI_COMM_WORLD, ierr)
          endif
      endif
      call chrtime(time,date,chtime,chdate)
c
c========================= Process Analysis End ======================
c
      if( iproc_id .EQ. 0 ) then
         if( date .EQ. wrtdate .AND. ABS(time-wrttim) .LT. 0.01 ) then
c
c-----Regular model grid files
c
           write(*,'(a20,$)') 'wrtcon ......'
           write(iout,'(a20,$)') 'wrtcon ......'

             do igrd = 1,ngrid
               nlayav = nlay(igrd)
               num_dims = 4
               if( .NOT. l3davg(igrd) ) then
                  nlayav = 1
                  num_dims = 3
               endif
               if( .NOT. lcdfout) then
                   call wrtcon(0,time,date,iavg(igrd),igrd,ncol(igrd),
     &                     nrow(igrd),nlayav,navspc,avcnc(iptrav(igrd)))
               else
                   if( igrd .eq. 1 )  then
                       action = 'Writing average concentration file '//
     &                                                     'for master grid.'
                   else
                       write(action,'(2A,I2)') 'Writing average ',
     &                                  'concntration file for grid: ',igrd 
                   endif
                   call ncf_wrt_data_tstep(action,iavg(igrd),navspc)
                   loutspec = .TRUE.
                   call ncf_wrt_data_species(action,iavg(igrd),igrd,num_dims,
     &                ncol(igrd),nrow(igrd),nlay(igrd),nlayav,nlayav,navspc,
     &                spavnam,loutspec,avcnc(iptrav(igrd)),height(iptr3d(igrd)))
               endif
             enddo
c          endif

           if( ldry .OR. lwet ) then
             do igrd = 1,ngrid
               if( .NOT. lcdfout ) then
                   call wrtdep(time,date,igrd,idep(igrd),ncol(igrd),nrow(igrd),
     &                             3*ndepspc+2,nspec,vdep(iptr1lay(igrd)),
     &                                            depfld(iptrdp(igrd)))
               else
                   if( igrd .eq. 1 )  then
                       action = 'Writing deposition file for master grid.'
                   else
                       write(action,'(2A,I2)') 'Writing deposition file ',
     &                                                  'for grid: ',igrd 
                   endif
                   ndps = 4*ndepspc
                   if( lixemis ) ndps = ndps + 2
                   call ncf_wrt_data_tstep(action,idep(igrd),ndps)
                   call ncf_wrt_data_dep(action,idep(igrd),igrd,ncol(igrd), 
     &                  nrow(igrd),nlay(igrd),3*ndepspc+2,nspec,depsp,
     &                             vdep(iptr1lay(igrd)),depfld(iptrdp(igrd)),
     &                              height(iptr3d(igrd)))
               endif
               if (lsrfmod) then
                 if( .NOT. lcdfout ) then
                    call wrtsrf(.TRUE.,igrd,time,date,ismout(igrd),
     &                       ncol(igrd),nrow(igrd),nsmspc,smspc,
     &                       solmas(iptrsm(igrd)),vegmas(iptrsm(igrd)))
                  else
                    if( igrd .eq. 1 )  then
                         action = 'Writing surface model file for master grid.'
                    else
                         write(action,'(2A,I2)') 'Writing surface model file',
     &                                                  'for grid: ',igrd 
                    endif
                    call ncf_wrt_data_tstep(action,ismout(igrd),nsmspc*2)
                    call ncf_wrt_data_srfmod(action,ismout(igrd),igrd,
     &                  ncol(igrd),nrow(igrd),nlay(igrd),nsmspc,
     &                    smspc,solmas(iptrsm(igrd)),vegmas(iptrsm(igrd)),
     &                            height(iptr3d(igrd)))
                  endif
               endif
             enddo
c
c======================== Source Apportion Begin =======================
c
             if( lptdepout ) then
                do igrd = 1,ngrid
                  if( .NOT. lcdfout ) then
                      call wrtdepsa(time,date,igrd,iowptdep(igrd),ncol(igrd),
     &                     nrow(igrd),notimespc,ptdryfld(ipsadep(igrd)),
     &                                         ptwetfld(ipsadep(igrd)))
                  else
                     if( igrd .eq. 1 )  then
                         action = 'Writing SA deposition file for master grid.'
                     else
                         write(action,'(2A,I2)') 'Writing SA deposition file',
     &                                                  ' for grid: ',igrd 
                     endif
                     call ncf_wrt_data_tstep(action,iowptdep(igrd),ntotspout)
                     call ncf_wrt_data_species(action,iowptdep(igrd),igrd,3,ncol(igrd),
     &                      nrow(igrd),nlay(igrd),1,1,ntotsp,ptdepname,loutsa,
     &                      ptdryfld(ipsadep(igrd)),height(iptr3d(igrd)))
                     call ncf_wrt_data_species(action,iowptdep(igrd),igrd,3,ncol(igrd),
     &                    nrow(igrd),nlay(igrd),1,1,ntotsp,ptdepname(ntotsp+1),
     &                    loutsa,ptwetfld(ipsadep(igrd)),height(iptr3d(igrd)))
                  endif
                enddo
             endif
             if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
                do igrd = 1,ngrid
                   if( .NOT. lcdfout ) then
                     call wrtsrf(lsrfmodrt,igrd,time,date,iowrtsrf(igrd),
     &                                  ncol(igrd),nrow(igrd),ntotsp,ptname,
     &                                  rtsolmas(ipsa2d(igrd)),rtvegmas(ipsa2d(igrd)) )
                   else
                     call ncf_wrt_data_tstep(action,iowrtsrf(igrd),ntotsp*2)
                     call ncf_wrt_data_species(action,iowrtsrf(igrd),igrd,3,ncol(igrd),
     &                      nrow(igrd),nlay(igrd),1,1,ntotsp,ptdepname,loutsa,
     &                      rtsolmas(ipsa2d(igrd)),height(iptr3d(igrd)))
                     call ncf_wrt_data_species(action,iowrtsrf(igrd),igrd,3,ncol(igrd),
     &                      nrow(igrd),nlay(igrd),1,1,ntotsp,ptdepname(ntotsp+1),loutsa,
     &                      rtvegmas(ipsa2d(igrd)),height(iptr3d(igrd)))
                   endif
                   nodes = ncol(igrd)*nrow(igrd)*ntotsp
                   call zeros(rtsolmas(ipsa2d(igrd)),nodes)
                   call zeros(rtvegmas(ipsa2d(igrd)),nodes)
               enddo
             endif
c
c======================== Source Apportion End =======================
c
c
           endif
c
           do igrd = 1,ngrid
             nlayav = nlay(igrd)
             if( .NOT. l3davg(igrd) ) nlayav = 1
             nodes = ncol(igrd)*nrow(igrd)*nlayav*navspc
             call zeros(avcnc(iptrav(igrd)),nodes)
             nodes = ncol(igrd)*nrow(igrd)*(3*ndepspc + 2)
             call zeros(depfld(iptrdp(igrd)),nodes)
             if( lptdepout ) then
                nodes = ncol(igrd)*nrow(igrd)*notimespc
                call zeros(ptdryfld(ipsadep(igrd)),nodes)
                call zeros(ptwetfld(ipsadep(igrd)),nodes)
             endif
           enddo
c
c-----Write the PiG sampling grid output
c
           if( lsample ) then
               do igrd = 1,nsample
                  if( .NOT. lcdfout ) then
                      call wrtsmp(.FALSE.,time,date,isample(igrd),
     &                         ncolsmp(igrd),nrowsmp(igrd),1,
     &                         navspc,smpcnc(ipsmp(igrd)))
                   else
                      write(action,'(2A,I2)') 'Writing data to sampling ',
     &                                    'grid file for sampling grid: ',igrd
                      call ncf_wrt_data_tstep(action,isample(igrd),navspc)
                      loutspec = .TRUE.
                      call ncf_wrt_data_species(action,isample(igrd),igrd,3,
     &                     ncolsmp(igrd),nrowsmp(igrd),nlay(igrd),
     &                          1,1,navspc,spavnam,loutspec,smpcnc(ipsmp(igrd)),
     &                                 height(iptr3d(igrd)))
                   endif
               enddo
c            endif
             do igrd = 1,nsample
               nodes = ncolsmp(igrd)*nrowsmp(igrd)*navspc
               call zeros(smpcnc(ipsmp(igrd)),nodes)
             enddo
           endif
c
c======================== Source Apportion Begin =======================
c
           if( ltrace .OR. lddm .OR. lhddm ) then
               if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) then
                  call wrrcprt(date,time)
               else
c
c   --- call routine to get the averages at the receptors and
c       write the receptor average file ----
c
                  do igrd=1,ngrid
                     nlayav = nlay(igrd)
                     if( .NOT. lsa_3davrg ) nlayav = 1
                     if( .NOT.( ltrace .OR. lddmcalc(igrd) ) ) cycle
                     call addrcp(igrd,ncol(igrd),nrow(igrd),nlayav,
     &                                   ntotsp,ptavrg(ipsa2d_avrg(igrd)))
                  enddo
                  if( ltrace ) call avgrcp(date,time)
                  if( lddm .OR. lhddm ) call avgrcpddm(date,time)
               endif
c
c
c   --- call routine to write the tracer gridded surface concentrations ---
c
              do igrd=1,ngrid
                if( .NOT.( ltrace .OR. lddmcalc(igrd) ) ) cycle
                nlayav = nlay(igrd)
                num_dims = 4
                if( .NOT. lsa_3davrg ) then
                   nlayav = 1
                   num_dims = 3
                endif
                if( .NOT. lcdfout ) then
                   call wsfcsa(igrd,date,time,ncol(igrd),nrow(igrd),nlayav,
     &                                        ntotsp,ptavrg(ipsa2d_avrg(igrd)))
                else
                   if( igrd .eq. 1 )  then
                       action = 'Writing tracer average concentration file '//
     &                                                     'for master grid.'
                   else
                       write(action,'(2A,I2)') 'Writing average ',
     &                                  'tracer concntration file for grid: ',igrd 
                   endif
                   call ncf_wrt_data_tstep(action,iowsfc(igrd),ntotspout)
                   call ncf_wrt_data_species(action,iowsfc(igrd),igrd,num_dims,ncol(igrd),
     &                      nrow(igrd),nlay(igrd),nlayav,nlayav,ntotsp,ptnameout,loutsa,
     &                         ptavrg(ipsa2d_avrg(igrd)),height(iptr3d(igrd)))
                endif
              enddo
c
c   --- call routine to re-initialize the running averages ---
c
              do igrd=1,ngrid
                 nlayav = nlay(igrd)
                 if( .NOT. lsa_3davrg ) nlayav = 1
                 if( .NOT.( ltrace .OR. lddmcalc(igrd) ) ) cycle
                 nodes=ncol(igrd)*nrow(igrd)*nlayav*ntotsp
                 call zeros(ptavrg(ipsa2d_avrg(igrd)),nodes)
              enddo
c
c   --- call routine to write the RTRAC/PiG sampling grid output
c
              if ((tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) .AND.
     &             lsample .AND. lsmptrc ) then
                   do igrd = 1,nsample
                       if( .NOT. lcdfout ) then
                           call wrtsmp(.true.,time,date,iowsmp(igrd),
     &                          ncolsmp(igrd),nrowsmp(igrd),1,nrtrac,
     &                          rtsmpcnc(iprtsmp(igrd)))
                       else
                           write(action,'(2A,I2)') 'Writing data to RTRAC ',
     &                            'sampling grid file for sampling grid: ',igrd
                           call ncf_wrt_data_tstep(action,iowsmp(igrd),nrtrac)
                           loutspec = .TRUE.
                           call ncf_wrt_data_species(action,iowsmp(igrd),igrd,3,
     &                          ncolsmp(igrd),nrowsmp(igrd),nlay(igrd),1,1,
     &                          nrtrac,ptname,loutspec,rtsmpcnc(iprtsmp(igrd)),
     &                          height(iptr3d(igrd)))
                       endif
                    enddo
c                endif
                 do igrd = 1,nsample
                   nodes = ncolsmp(igrd)*nrowsmp(igrd)*nrtrac
                   call zeros(rtsmpcnc(iprtsmp(igrd)),nodes)
                 enddo
              endif
           endif
c
c========================= Source Apportion End ========================
c
c========================= Process Analysis Begin ======================
c
c-----Get final concentration
c
           if( lipr .OR. lirr ) then
              if( lipr .AND. .NOT. lmpi ) then
                linit = .FALSE.
                do igrd = 1,ngrid
                   call initipr(linit,iproc_id,igrd,nspec,ncol(igrd),
     &                        nrow(igrd), nlay(igrd),conc(iptr4d(igrd)))
                enddo
              endif
c
c-----Write PA results
c
              if( lipr ) call wrtipr(date,time)
              if( lirr ) call wrtirr(date,time)
c
c   --- call routine to zero out all Process Analysis data structures ---
c
              if( .NOT. lmpi ) call pazero( )
c
c-----Get initial concentration for next loop
c
              if( lipr .AND. .NOT. lmpi ) then
                 linit = .TRUE.
                 do igrd = 1,ngrid
                    call initipr(linit,iproc_id,igrd,nspec,ncol(igrd),
     &                        nrow(igrd),nlay(igrd),conc(iptr4d(igrd)))
                 enddo
              endif
           endif
c
c========================= Process Analysis End ========================
c
           write(*,'(a)') '   Done'
           write(iout,'(a)') '   Done'
           call flush(6)
           call flush(iout)
c
c-----Write PiG restart file and diagnostics
c
           if( ipigflg .NE. 0 ) then
             write(*,'(a20,$)') 'wrtpig ......'
             write(iout,'(a20,$)') 'wrtpig ......'
             call wrtpig(date,time,begdate,begtim)
c            call pigdiag(idiag,chtime,chdate,1,'                    ')
             call pigmscl(nspec,ngrid,chtime,chdate,idiag,
     &                                                  pigdump,pgmserr)
             write(*,'(a)') '   Done'
             write(iout,'(a)') '   Done'
             call flush(6)
             call flush(iout)
           endif
c
c-----Write model mass
c
           do igrd = 1,ngrid
             call wrtmass(igrd,chdate,chtime,1)
           enddo
c
c-----Flush file units
c
           call flush(iout)
           call flush(idiag)
           call flush(imass)
           ncf_cur_tstep = ncf_cur_tstep + 1
         endif
      endif
c
      if( lmpi .AND. iproc_id .GT. 0 ) then
         if( date .EQ. wrtdate .AND. ABS(time-wrttim) .LT. 0.01 ) then
             do igrd = 1,ngrid
               call newgrid(igrd)
               nlayav = nlay(igrd)
               if( .NOT. l3davg(igrd) ) nlayav = 1
               nodes = mxp*myp*nlayav*navspc
               call zeros(avcnc(iptrav(igrd)),nodes)
               nodes = mxp*myp*(3*ndepspc + 2)
               call zeros(depfld(iptrdp(igrd)),nodes)
               call nodemass(igrd)
c
c========================= Source Apportion Begin ========================
c
               if( ltrace .OR. lddm .OR. lhddm ) then
                 nlayav = nlay(igrd)
                 if( .NOT. lsa_3davrg ) nlayav = 1
                 nodes=mxp*myp*nlayav*ntotsp
                 call zeros(ptavrg(ipsa2d_avrg(igrd)),nodes)
                 if( lptdepout ) then
                    nodes = mxp*myp*notimespc
                    call zeros(ptdryfld(ipsadep(igrd)),nodes)
                    call zeros(ptwetfld(ipsadep(igrd)),nodes)
                 endif
               endif
c
               if( (tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) .AND.
     &             .not.lsrfmodrt ) then
                  nodes = mxp*myp*ntotsp
                  if (lsrfmod) then
                     call zeros(rtsolmas(ipsa2d(igrd)),nodes)
                     call zeros(rtvegmas(ipsa2d(igrd)),nodes)
                  endif
               endif
c
c========================= Source Apportion End ========================
c
             enddo 
             if( ipigflg .NE. 0 ) then
                 call zeros(nkill,9)
                 call zeros(nage,ngrid)
                 call zeros(pigage,ngrid)
             endif
             if( ipigflg .NE. 0 .AND. lsample ) then
                 do ismp=1,nsample
                   nodes = ncolsmp(ismp)*nrowsmp(ismp)*navspc
                   call zeros(smpcnc(ipsmp(ismp)),nodes)
                   if( (tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC)
     &                                           .AND. lsmptrc ) then
                       nodes = ncolsmp(ismp)*nrowsmp(ismp)*nrtrac
                       call zeros(rtsmpcnc(iprtsmp(ismp)),nodes)
                   endif
                 enddo
             endif
         endif
      endif
c
      if( date .EQ. wrtdate .AND. ABS(time-wrttim) .LT. 0.01 ) then
         whr = aint(wrttim/100.)
         wmn = amod(wrttim,100.)
         wrttim = 100.*(whr + aint((wmn + dtout)/60.)) + 
     &                                     amod((wmn + dtout),60.)
         if (wrttim.ge.2400.) then
          wrttim = wrttim - 2400.
          wrtdate = wrtdate + 1
          if( MOD(wrtdate,1000) .GT. 365 ) then
              if( MOD(INT(wrtdate/1000),4) .EQ. 0 ) then
                 if( MOD(wrtdate,1000) .EQ. 367 )
     &                     wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
              else
                 wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
              endif
           endif
         endif
      endif
c
c-----Check for end of simulation
c
      if (date.lt.enddate) goto 100
      if (date.eq.enddate .and. time.lt.endtim - 0.01) goto 100
c
c------------------  End main time-integration loop  -------------------
c
c-----Send the instantaneous fields back to master ---
c
      if( lmpi ) then
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         call conc_update(numprocs, iproc_id)
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
      endif
c
c-----Write final instantaneous restart files
c
      if( iproc_id .EQ. 0 ) then
          call wrtcon(1,time,date,iconc,1,ncol(1),nrow(1),nlay(1),
     &                                                   nspec,conc(1))
          if( ngrid .GT. 1) call wrfgcon(date,time)
          if( lcdfout ) then
             do i=1,ngrid
                if( igrd .eq. 1 )  then
                  action = 'Closing average concentration file '//
     &                                                     'for master grid.'
                else
                  write(action,'(2A,I2)') 'Closing average ',
     &                                  'concntration file for grid: ',i
                endif
                if( ldry .OR. lwet ) then 
                   if( igrd .eq. 1 )  then
                     action = 'Closing deposition file for master grid.'
                   else
                     write(action,'(2A,I2)') 'Closing deposition ',
     &                                            'file for grid: ',i
                   endif
                endif
             enddo
         endif
c
c======================== Source Apportion Begin =======================
c
          if( ltrace .OR. lddm .OR. lhddm ) then
             if( lddmcalc(1) ) call wconsa(date,time,ncol(1),nrow(1),
     &                                        nlay(1),ntotsp,ptconc(1))
             if( ngrid .GT. 1 ) call wfconsa(date,time)
          endif
      endif
c
c========================= Source Apportion End ========================
c
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      write(iout,'(/,a,a8,1x,a8,/)') 'Date/time: ',chdate,chtime
      write(iout,'(a)')'END SIMULATION'
      write(iout,'(a,i10)') 'TOTAL MASTER GRID TIME STEPS: ', nsteps
      if( iproc_id .EQ. 0 ) then
         write(*,'(/,a,a8,1x,a8,/)')'Date/time: ',chdate,chtime
         write(*,'(a)')'END SIMULATION'
         write(*,'(a,i10)') 'TOTAL MASTER GRID TIME STEPS: ', nsteps
      end if
c
c  --- call routine to close any open NetCDF datasets ---
c
      if( lcdfout ) call ncf_closefiles()
c
      if (lmpi) call MPI_Finalize(ierr)
      stop
      end

