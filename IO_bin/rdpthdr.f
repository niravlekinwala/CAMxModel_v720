      subroutine rdpthdr(begtim,begdate,endtim,enddate)
      use filunit
      use grid
      use chmstry
      use bndary
      use ptemiss
      use tracer
c 
c----CAMx v7.20 220430
c 
c     RDPTHDR reads the header of binary point source emissions file,
c     initializes time-invariant point source variables, and maps the point
c     source species list to the internal CAMx species list
c                           
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        10/24/01  Removed BSWAP and converted integer strings to character*4
c        11/06/01  Input dates are now Julian
c        02/09/02  Added code to handle end of year dates
c        12/21/13  Added compact point source file
c        12/07/14  Revised for VBS emissions
c        01/08/16  Updated for Revised VBS emissions
c        11/22/17  Added checks for PFE/PMN/PK/PCA/PMG emissions (disabled)
c        03/15/21  Added check for number of species in file
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
      integer ifle
c
      integer, parameter :: num_elem = 5
      integer :: idx_elem(num_elem)
      integer, allocatable, dimension(:) :: idxpnt
      integer idx_start, n
      logical, allocatable, dimension(:) :: lfound

c
      data ptfil /'PTSOURCE  '/
      data ptcompact /'PTSOURCECP'/
c
c-----Entry point
c
c-----Initialize the mapping array ---
c
      do ifle=1,npoint_files
c
c-----skip if this is a NetCDF file ---
c
        if( is_netcdf_iptem(ifle) ) cycle
c
c-----Read 1st PT header record and check inputs 
c             
        rewind(iptem(ifle))
        read(iptem(ifle)) ifile,note,nseg,nptspc(ifle),idat1,tim1,idat2,tim2
        if( nptspc(ifle) .GT. MXSPEC ) then
             write(iout,'(//,a)') 'ERROR in RDPTHDR:'
             write(iout,*) 'Number of point source emissions species exceeds max.'
             write(iout,*) 'Increase the parameter MXSPEC to at least: ',nptspc(ifle)
             call camxerr()
        endif
        allocate( idxpnt(nptspc(ifle)) )
        allocate( lfound(nptspc(ifle)) )
        idxpnt = 0
c
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
          write(iout,'(//,a)') 'ERROR in RDPTHDR:'
          write(iout,*)'PT input file is not labelled PTSOURCE or PTSOURCECP' 
          call camxerr()
        endif   
        tim1 = 100.*tim1
        tim2 = 100.*tim2
        if (idat1.gt.begdate) then 
          write(iout,'(//,a)') 'WARNING in RDPTHDR:'
          write(iout,*)'PT start date > simulation start date' 
          write(iout,*)'  PT file: ',idat1 
          write(iout,*)'Sim start: ',begdate 
          if (.not.le1day) then
            write(iout,*)'CAMx expecting day-specific emissions: Stopping'
            call camxerr()
          endif
        elseif (idat1.eq.begdate .and. tim1.gt.begtim) then 
          write(iout,'(//,a)') 'WARNING in RDPTHDR:'
          write(iout,*)'PT start time > simulation start time' 
          write(iout,*)'  PT file: ',tim1 
          write(iout,*)'Sim start: ',begtim 
          call camxerr()
        elseif (idat2.lt.enddate) then 
          write(iout,'(//,a)') 'WARNING in RDPTHDR:'
          write(iout,*)'PT end date < simulation end date' 
          write(iout,*)'PT file: ',idat2
          write(iout,*)'Sim end: ',enddate 
          if (.not.le1day) then
            write(iout,*)'CAMx expecting day-specific emissions: Stopping'
            call camxerr()
          endif
        elseif (idat2.eq.enddate .and. tim2.lt.endtim) then 
          write(iout,'(//,a)') 'ERROR in RDPTHDR:'
          write(iout,*)'PT end time < simulation end time' 
          write(iout,*)'PT file: ',tim2
          write(iout,*)'Sim end: ',endtim 
          call camxerr()
        endif 
c 
c-----Read 2nd PT header
c 
        read(iptem(ifle)) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz 
c 
c-----Read 3rd & 4th PT header 
c 
        read(iptem(ifle))
        read(iptem(ifle)) ((ptspec(n,l),n=1,10),l=1,nptspc(ifle)) 
c
        do 15 lpt = 1,nptspc(ifle) 
           write(ptspc,'(10a1)') (ptspec(n,lpt),n=1,10) 
           if (ptspc.eq.'HNO2      ') ptspc = 'HONO      '
           if (ptspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        ptspc = 'FORM      '
c
c  --- find species in list ---
c
           lfound(lpt) = .FALSE.
           do ispc=1,nspec
              if(ptspc .EQ. spname(ispc)) lfound(lpt) = .TRUE.
           enddo
           if( lvbs .AND. LVBSPREPROC ) then
              if( ptspc .eq. 'IVOA      ' ) lfound(lpt) = .TRUE.
              if( ptspc .eq. 'IVOB      ' ) lfound(lpt) = .TRUE.
              if( ptspc .eq. 'POA_OP    ' ) lfound(lpt) = .TRUE.
              if( ptspc .eq. 'POA_BB    ' ) lfound(lpt) = .TRUE.
           endif
           if( .NOT. lfound(lpt) ) cycle
           do k=1,nemspc
             if(ptspc .EQ. emspcname(k) ) idxpnt_files(ifle,lpt) = k
           enddo
           if( idxpnt_files(ifle,lpt) .EQ. 0 ) then
              nemspc = nemspc + 1
              emspcname(nemspc) = ptspc
              idxpnt_files(ifle,lpt) = nemspc
           endif
c
c  --- check if species is not modelled ---
c   
           do 20 l = 1,nspec 
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
               lemmap(l) = idxpnt_files(ifle,lpt)
               idxpnt(lpt) = l
               write(idiag,'(2(a,i5,2x,a))')
     &                   'Point source species ',lpt,ptspc, 
     &                   ' mapped to model species ',l,spname(l) 
             endif 
 20        continue
           if( idxpnt(lpt) .LE. 0 ) write(idiag,*)'PT species: ',ptspc,' not modeled'
 15      continue 
c
c-----Read time invariant data (assume nseg=1)
c     check number of sources against max 
c
         read(iptem(ifle)) idum,nptsrc_files(ifle)
         if( allocated(idxpnt) ) deallocate( idxpnt )
         if( allocated(lfound) ) deallocate( lfound )
      enddo
c
c----Call routine to allocate the arrays 
c
      call alloc_ptemiss(nspec,ngrid,nptsrc)
c
c  --- read rest of each file ---
c
      do ifle=1,npoint_files
c
c-----skip if this is a NerCDF file ---
c
         if( is_netcdf_iptem(ifle) ) cycle
         idx_start = idx_start_pts(ifle)
         read(iptem(ifle)) (xloc(n),yloc(n),hstk(n),dstk(n),tstk(n),
     &                     vstk(n),n=idx_start+1,idx_start+nptsrc_files(ifle))
c
         do 30 n = idx_start+1,idx_start+nptsrc_files(ifle)
           lpiglet(n) = .false.
           if (dstk(n).lt.0.) lpiglet(n) = .true.
           vstk(n) = vstk(n)/3600.
c
c========================= Source Apportion End ========================
c
  30     continue
      enddo
      write(idiag,*)
c
      return
      end
