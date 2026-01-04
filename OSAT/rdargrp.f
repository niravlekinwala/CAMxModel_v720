C***** RDARGRP
c
      subroutine rdargrp(igrid,ndate,ttime,nox,noy,nlay_ems,
     &                                          numcls,igroup,emscls)
      use filunit
      use grid
      use bndary
      use camxcom
      use tracer
      implicit none
c
c      Copyright 1996 - 2022
c     Ramboll
c
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c      10/24/01  Removed BSWAP and converted integer strings to character*4
c      05/01/03  Time span of emissions must now match emiss update interval
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
      integer igrid
      integer ndate
      real    ttime
      integer nox
      integer noy
      integer nlay_ems
      integer numcls
      integer igroup
      real    emscls(numcls,nox,noy,nlay_ems)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  cname
      character*4   iname(10)
      integer       ispc, num_emsfiles, idxfile, itmp, jtmp
      integer       ibgdat, iendat, iounit, idx, iseg, i, j
      integer       ichkems, ichktm1, ichktm2, icls, buffer_offset
      real          bgtim, edtim
      logical       lfound, use_this_file
c
      real,allocatable, dimension(:,:) :: emsgrd
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- initalize the array ---
c
      emscls = 0.
c
c   --- allocate loacal array ---
c
      allocate( emsgrd(nox,noy) )
c
c   --- skip if filename not supplied ---
c
      if( .NOT. larsrc ) goto 9999
      num_emsfiles = num_iortem(igrid,igroup)
      if( igroup .EQ. 0 ) num_emsfiles = nemiss_files(igrid)
c
c   --- loop over all files ---
c
      do 10 idxfile=1,num_emsfiles
        use_this_file = .TRUE.
c
c   --- set the unit number for surface emissions file ---
c
        if( igroup .EQ. 0 ) then
            iounit = iarem(igrid,idxfile)
            buffer_offset = buffer_offset_iarem(igrid,idxfile)
            if( is_netcdf_iarem(igrid,idxfile) ) then
                use_this_file = .FALSE.
            else
                write(fname,'(A,I3)') 'EMISSIONS -- UNIT ',iarem(igrid,idxfile)
c
c    --- if emissions is regular model emissions file, backup one hour ---
c
                do i=1,nspcem(igrid,igroup,idxfile)
                    backspace(iounit)
                enddo
                backspace(iounit)
            endif
        else
            iounit = iortem(igrid,igroup,idxfile)
            fname = temfil(igrid,igroup,idxfile)
            if( is_netcdf_iortem(igrid,igroup,idxfile) ) use_this_file = .FALSE.
            buffer_offset = buffer_offset_iortem(igrid,igroup,idxfile)
        endif
        if( .NOT. use_this_file ) cycle
c
c   --- read the date and time, again ---
c
        lfound = .FALSE.
  111   continue
        read(iounit) ibgdat, bgtim, iendat, edtim
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
           write(iout,'(//,a)')'ERROR in RDARGRP:'
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
c   --- read the emissions for this hour ---
c
        do 20 ispc=1,nspcem(igrid,igroup,idxfile)
            read(iounit,ERR=7000) iseg, (iname(i),i=1,10), 
     &                         ((emsgrd(i,j),i=1+buffer_offset,nox-buffer_offset),
     &                                            j=1+buffer_offset,noy-buffer_offset)
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
c   --- if the species is a not modelled or not a VOC species skip it ---
c
            idx = idxems(igrid,igroup,idxfile,ispc)
            if( idx .LE. 0 ) goto 20
            if( .NOT. lusespc(idx) ) goto 20
c
c   --- load the species into the cells in tracer array ----
c
            do j=1+buffer_offset,noy-buffer_offset
               jtmp = j-buffer_offset
               do i=1+buffer_offset,nox-buffer_offset
                  do icls=1,ntrcls
                     itmp = i-buffer_offset
                     emscls(icls,i,j,1) = emscls(icls,i,j,1) +
     &                      (emsgrd(itmp,jtmp)/(60.*dtems)) * trspmap(idx,icls)
                  enddo
               enddo
           enddo
c
c   --- next species ---
c
  20    continue
c
c   --- if the correct hour has not been found, 
c       go back and read some more else read next file ---
c
        if( .NOT. lfound ) then
           goto 111
        else
           if( idxfile .LT. num_emsfiles ) then
              goto 10
           else
              goto 9999
           endif
        endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222   continue
        if( le1day ) then
            rewind(iounit)
            read(iounit)
            read(iounit)
            read(iounit)
            read(iounit)
            goto 111
        else
           goto 7001
        endif
c
c  --- next file in grid/group ---
c
      deallocate( emsgrd )
 10   continue
c
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDARGRP:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'Reading emissions ',
     &   'after hour ',ibgdat, bgtim,' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDARGRP:'
      write(iout,'(/,1X,2A,I8.5,F4.1,2A)') 'Premature end-of-file',
     &              ' in emissions file after hour ',ibgdat, bgtim,
     &              ' in file: ',fname(:istrln(fname))
      call camxerr()
c
 9999 continue
c
      return
      end
