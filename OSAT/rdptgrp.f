C***** RDPTGRP
c
      subroutine rdptgrp(ndate,ttime,igroup,numcls,emscls,izcel)
      use filunit
      use grid
      use camxcom
      use chmstry
      use ptemiss
      use tracer
c
c
c      Copyright 1996 - 2022
c     Ramboll
c
c
c----CAMx v7.20 220430
c
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     08/18/99  --gwilson--   Added code to implement the override
c                             flag for the source area of point sources
c     10/24/01  Removed BSWAP and converted integer strings to character*4
c     07/05/02  Changed to account for new type of the PiG flag
c     05/01/03  Time span of emissions must now match emiss update interval
c     09/25/03  Significant changes to handle PSAT
c     03/01/16  Fixed bug in point source override and compact ptsrce file
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer ndate
      real    ttime
      integer igroup
      integer numcls
      real    emscls(MXTRCLS,*)
      integer izcel(*)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, idx, iseg, idxfile
      integer       npoint, ispc, i, idum, idxcp, idcompact(MXPTSRC)
      integer       num_ptsfiles, idx_start
      real          bgtim, edtim, rdum
      logical       lfound, lfirst_pass
c
      integer, allocatable, dimension(:) :: izcel_in
      real,    allocatable, dimension(:) :: emspts_in
      real,    allocatable, dimension(:) :: emspts
      real,    allocatable, dimension(:) :: flowrat_in
      real,    allocatable, dimension(:) :: effph_in
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- allocate local arrays ----
c
      allocate( izcel_in(nptsrc) )
      allocate( emspts_in(nptsrc) )
      allocate( emspts(nptsrc) )
      allocate( flowrat_in(nptsrc) )
      allocate( effph_in(nptsrc) )
c
c   --- skip if filename not supplied ---
c
      num_ptsfiles = num_iortpt(igroup)
      if( igroup .EQ. 0 ) num_ptsfiles = npoint_files
      do idxfile=1,num_ptsfiles
         if( .NOT. ltptfl(igroup,idxfile) ) cycle
         idx_start = idx_start_sapts(igroup,idxfile)
c
c   --- set the unit number for file ---
c
         if( igroup .EQ. 0 ) then
             if( is_netcdf_iptem(idxfile) ) cycle
             iounit = iptem(idxfile)
             write(fname,'(A,I3)') 'PTSOURCE -- UNIT ',iptem(idxfile)
             do i=1,nspcpt(igroup,idxfile)
                 backspace(iounit)
             enddo
             backspace(iounit)
             backspace(iounit)
             backspace(iounit)
         else
             if( is_netcdf_iortpt(igroup,idxfile) ) cycle
             iounit = iortpt(igroup,idxfile)
             fname = tptfil(igroup,idxfile)
         endif
         lfirst_pass = .TRUE.
c
c  --- gary wants this flush ---
c
         call flush(6)
c
c   --- read the date and time, again ---
c
         lfound = .FALSE.
  111    continue
         read(iounit,END=222) ibgdat, bgtim, iendat, edtim
         if(NINT(edtim) .EQ. 0) then
           edtim = 24.
           iendat = ibgdat
         endif
         ichktm1 = NINT( 1000*(bgtim) )
         if( le1day ) then
            ichktm2 = NINT( 1000*(edtim) )
         else
            ichktm2 = NINT( 1000*(edtim)+24000*(iendat-ibgdat) )
         endif
         if( NINT(edtim) .EQ. 0 ) ichktm2 = 24000
         ichkems = NINT( 1000*(dtems/60.) )
         if( (ichktm2 - ichktm1) .NE. ichkems ) then
             write(iout,'(//,a)')'ERROR in RDPTGRP:'
             write(iout,*) 'Time interval in surface emissions file does'
             write(iout,*)  ' not match emissions update time interval.'
             write(iout,*) '   Beginning Date/Time (Hour): ',ibgdat,bgtim
             write(iout,*) '   Ending Date/Time    (Hour): ',iendat,edtim
             write(iout,*) '   Emiss Input interval(Hour): ',dtems/60.
             call camxerr()
         endif
         bgtim = 100.*aint(bgtim) + 60.*amod(bgtim,1.)
         edtim = 100.*aint(edtim) + 60.*amod(edtim,1.)
c
c   --- read the number of points and point locations ---
c
          read(iounit,ERR=7000) iseg, npoint
          if( npoint .LE. 0 ) cycle
          read(iounit,ERR=7000,END=7001)  (idcompact(i), idum,
     &             izcel_in(i), flowrat_in(i), effph_in(i),i=1,npoint)
          do i=1,npoint
            idxcp = i
            if( lcompactpt(igroup,idxfile) ) idxcp = idcompact(i)
            izcel(i+idx_start) = izcel_in(i)
            flowrat(i+idx_start) = flowrat_in(i)
            effph(i+idx_start) = effph_in(i)
            if( flowrat(i) .LT. 0. .AND. effph(i) .LT. 0. .AND.
     &                    ABS(effph(i)) .LT. ABS(flowrat(i)) ) goto 7003
          enddo
c
c   --- read the emissions for this hour ---
c
          do 10 ispc=1,nspcpt(igroup,idxfile)
             read(iounit,ERR=7000,END=7001) iseg, (iname(i),i=1,10),
     &                                         (emspts_in(i),i=1,npoint)
c
c   --- if date and time does not match this hour, skip this record ---
c
              if( le1day ) then
                  if( bgtim .NE. ttime ) goto 10
              else
                  if( ndate .NE. ibgdat .OR. bgtim .NE. ttime ) goto 10
              endif
              lfound  = .TRUE.
c
c   --- if the species is a not modelled or not a tracer species skip it ---
c
              idx = idxpts(igroup,idxfile,ispc)
              if( idx .LE. 0 ) goto 10
c
c   --- put emissions into regular model array ---
c
              do i=1,npoint
                idxcp = i
                if( lcompactpt(igroup,idxfile) ) idxcp = idcompact(i)
                ptemis(idxcp+idx_start,lemmap(idx)) = emspts_in(i)/(60.*dtems)
              enddo
c
c  --- skip the rest if not an SA species ---
c
              if( .NOT. lusespc(idx) ) goto 10
c
c   --- convert to PPM (PPMC) ----
c
              do i=1,npoint
                idxcp = i
                if( lcompactpt(igroup,idxfile) ) idxcp = idcompact(i)
                emspts(idxcp) = emspts_in(i)/(60.*dtems)
c
c   --- put into regular model array ---
c
                ptemis(i+idx_start,lemmap(idx)) = emspts_in(i)/(60.*dtems)
c
c   --- load into the tracer emissions array ---
c
                 do 20 icls=1,ntrcls
c
c   --- if this is the NOx emissions tracer and this is a 
c       PiG source, skip it (PiG treated elsewhere) ---
c
c   --- gary wants this call to flush ---
c
                     call flush(6)
                     if( lpigsa(idxcp+idx_start) ) goto 20
                     emscls(icls,idxcp+idx_start) = 
     &                   emscls(icls,idxcp+idx_start) + 
     &                          emspts(idxcp) * trspmap(idx,icls)
   20            continue
              enddo
c
c   --- next species ---
c
  10      continue
c
c   --- if the correct hour has not been found,
c       go back and read some more else read next file ---
c
          if( .NOT. lfound ) then
             goto 111
          else
             cycle
          endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222    continue
         if( le1day .AND. lfirst_pass ) then
             lfirst_pass = .FALSE.
             rewind(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             goto 111
         else
             goto 7001
         endif
         cycle
      enddo
      goto 9999
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDPTGRP:'
      write(iout,'(/,1X,A,I8.5,F8.1,2A)') 
     &      'Reading emissions after hour ',ibgdat, bgtim,' in file: ',
     &                                            fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDPTGRP:'
      write(iout,'(/,1X,3A)') 'Premature end-of-file reading ',
     &                  'emissions from file: ',fname(:istrln(fname))
      write(iout,'(1X,3A)') 'Make sure the file congtains all of the',
     &                  ' hours in the simulation.'
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDPTGRP:'
      write(iout,'(/,1X,A,I10,2A)') 'Number of points: ',npoint,
     &     ' is not consistent with regular emissions in file: ',
     &                                            fname(:istrln(fname))
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)')'ERROR in RDPTGRP:'
      write(iout,*) 'Invalid values found for flow rate and plume height.'
      write(iout,*) 'The plume distribution override is triggered. But the '
      write(iout,*) 'value for flow rate (base) is larger than value for ',
     &                      'plume height (top).'
      call camxerr()
c
 9999 continue
c
c   --- deallocate local arrays ----
c
      deallocate( izcel_in )
      deallocate( emspts_in )
      deallocate( emspts )
      deallocate( flowrat_in )
      deallocate( effph_in )
c
      return
      end
