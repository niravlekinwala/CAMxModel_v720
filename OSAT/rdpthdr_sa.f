      subroutine rdpthdr_sa(begtim,begdate,endtim,enddate)
      use filunit
      use grid
      use chmstry
      use bndary
      use ptemiss
      use tracer
c 
c----CAMx v7.20 220430
c 
c     RDPTHDR_SA reads the header of binary point source emissions file
c     in each soruce group, initializes time-invariant point source 
c     variables, and maps the point source species list to the internal 
c     CAMx species list
c                           
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c 
c     Input arguments: 
c        begtim              model start time (HHMM) 
c        begdate             model start date (YYJJJ) 
c        endtim              model end time (HHMM) 
c        enddate             model end date (YYJJJ)
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
      include 'vbs.inc'
c
      integer begdate
      integer enddate
c
      character*10 ptfil, infil, ptspc, ptcompact
      character*4  ifile(10), note(60)
c
      character*4 ptspec(10,MXSPEC)
c
      logical l_exist_ivoc, l_exist_poa_xx, lfound
      integer ifle, igroup, iounit
c
      integer, parameter :: num_elem = 5
      integer :: idx_elem(num_elem)
      integer idx_start, numpts, n
c
      data ptfil /'PTSOURCE  '/
      data ptcompact /'PTSOURCECP'/
c
c-----Entry point
c
c-----Loop over groups and files in each group ---
c
       numpts = 0
       do igroup=1,ngroup 
       do ifle=1,num_iortpt(igroup)
        iounit = iortpt(igroup,ifle)
c
c-----skip if this is a NetCDF file ---
c
        if( is_netcdf_iortpt(igroup,ifle) ) cycle
c
c-----Read 1st PT header record and check inputs 
c             
        rewind(iounit)
        read(iounit) ifile,note,nseg,nspcpt(igroup,ifle),idat1,tim1,idat2,tim2
        if(INT(tim2) .EQ. 24) then
          tim2 = 0.
          idat2 = idat2 + 1
          if( MOD(idat2,1000) .GT. 365 ) then
             if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
                if( MOD(idat2,1000) .EQ. 367 ) 
     &                       idat2 = (INT(idat2/1000)+1)*1000 + 1
             else
                idat2 = (INT(idat2/1000)+1)*1000 + 1
             endif
          endif
        endif
        write(infil,'(10a1)') (ifile(n),n=1,10) 
        if (infil.ne.ptfil .AND. infil .ne. ptcompact) then 
          write(iout,'(//,a)') 'ERROR in RDPTHDR_SA:'
          write(iout,*)'PT input file is not labelled PTSOURCE or PTSOURCECP' 
          call camxerr()
        endif   
        tim1 = 100.*tim1
        tim2 = 100.*tim2
        if (idat1.gt.begdate) then 
          write(iout,'(//,a)') 'WARNING in RDPTHDR_SA:'
          write(iout,*)'PT start date > simulation start date' 
          write(iout,*)'  PT file: ',idat1 
          write(iout,*)'Sim start: ',begdate 
          if (.not.le1day) then
            write(iout,*)'CAMx expecting day-specific emissions: Stopping'
            call camxerr()
          endif
        elseif (idat1.eq.begdate .and. tim1.gt.begtim) then 
          write(iout,'(//,a)') 'WARNING in RDPTHDR_SA:'
          write(iout,*)'PT start time > simulation start time' 
          write(iout,*)'  PT file: ',tim1 
          write(iout,*)'Sim start: ',begtim 
          call camxerr()
        elseif (idat2.lt.enddate) then 
          write(iout,'(//,a)') 'WARNING in RDPTHDR_SA:'
          write(iout,*)'PT end date < simulation end date' 
          write(iout,*)'PT file: ',idat2
          write(iout,*)'Sim end: ',enddate 
          if (.not.le1day) then
            write(iout,*)'CAMx expecting day-specific emissions: Stopping'
            call camxerr()
          endif
        elseif (idat2.eq.enddate .and. tim2.lt.endtim) then 
          write(iout,'(//,a)') 'ERROR in RDPTHDR_SA:'
          write(iout,*)'PT end time < simulation end time' 
          write(iout,*)'PT file: ',tim2
          write(iout,*)'Sim end: ',endtim 
          call camxerr()
        endif 
c 
c-----Read 2nd PT header
c 
        read(iounit) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz 
c 
c-----Read 3rd & 4th PT header 
c 
        read(iounit)
        read(iounit) ((ptspec(n,l),n=1,10),l=1,nspcpt(igroup,ifle)) 
c 
c-----Map PT species to model species 
c 
        do 15 lpt = 1,nspcpt(igroup,ifle) 
           write(ptspc,'(10a1)') (ptspec(n,lpt),n=1,10) 
           if (ptspc.eq.'HNO2      ') ptspc = 'HONO      '
           if (ptspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        ptspc = 'FORM      '
c
c  --- check if species is not modelled ---
   
           do l = 1,nspec 
             if( ptspc .EQ. spname(l) ) then 
               write(idiag,'(2(a,i5,2x,a))')
     &                   'Point source species ',lpt,ptspc, 
     &                   ' mapped to model species ',l,spname(l) 
               lfound = .FALSE.
               do k=1,nemspc
                 if( ptspc .EQ. emspcname(k) ) lfound = .TRUE.
               enddo
               if( .NOT. lfound ) then
                 nemspc = nemspc + 1
                 emspcname(nemspc) = ptspc
                 lemmap(l) = nemspc
               endif 
             endif
           enddo
   15   continue 
        if( lvbs .AND. LVBSPREPROC ) then
           if( ptspc .EQ. 'POA_OP    ') then
               lemmap(kpap_c(0)) = lpt
               write(idiag,'(a,i5,2x,2a)') 'Point source species ',
     &                   lpt,ptspc,' mapped to model species PAP0-PAP5 (VBS)'
           endif
           if( ptspc .EQ. 'POA_BB    ') then
               lemmap(kpfp_c(0)) = lpt ! temporarily assign to PFP0
               write(idiag,'(a,i5,2x,2a)') 'Point source species ',
     &               lpt,ptspc,' mapped to model species PFP0-PFP5 (VBS)'
           endif
        endif
        if( ptspc .EQ. spname(l) ) then
          write(idiag,'(2(a,i5,2x,a))')
     &                 'Point source species ',lpt,ptspc,
     &                 ' mapped to model species ',l,spname(l)
        endif
c
c-----Read time invariant data (assume nseg=1)
c     check number of sources against max 
c
         read(iounit) idum,nptsrc_safile(igroup,ifle)
         idx_start_sapts(igroup,ifle) = numpts
         numpts = numpts + nptsrc_safile(igroup,ifle)
c
c---- next file and group
c
      enddo
      enddo
c
c----Call routine to allocate the arrays 
c
      call alloc_ptemiss(nspec,ngrid,nptsrc)
c
c  --- read rest of each file ---
c
      do igroup=1,ngroup 
      do ifle=1,num_iortpt(igroup)
         iounit = iortpt(igroup,ifle)
c
c-----skip if this is a NerCDF file ---
c
         if( is_netcdf_iortpt(igroup,ifle) ) cycle
         idx_start = idx_start_sapts(igroup,ifle)
         read(iounit) (xloc(n),yloc(n),hstk(n),dstk(n),tstk(n),
     &                     vstk(n),n=idx_start+1,idx_start+nptsrc_safile(igroup,ifle))
         do 30 n = idx_start+1,idx_start+nptsrc_safile(igroup,ifle)
           lpiglet(n) = .false.
           if (dstk(n).lt.0.) then
               lpiglet(n) = .true.
           endif
           vstk(n) = vstk(n)/3600.
c
c========================= Source Apportion End ========================
c
  30     continue
      enddo
      enddo
      write(idiag,*)
c
      return
      end
