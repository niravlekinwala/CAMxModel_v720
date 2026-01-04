c***** RDARDDM.F
c
      subroutine rdarddm(igrid,ndate,ttime,nox,noy,nlay_ems,nddmspc,emisddm)
      use filunit
      use grid
      use chmstry
      use bndary
      use camxcom
      use tracer
c
c----CAMx v7.20 220430
c
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions for the DDM process
c   and fills the approproate arrays.  The emissions file for one grid
c   but each emissions groups is read.
c    Argument descriptions:
c     Outputs:
c       emisddm   R    array of emissions for DDM speices
c     Inputs:
c       igrid     I    grid number
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       nox       I    number of columns in grid
c       noy       I    number of rows in grid
c       nlay_ems  I    number of layers in emissions files
c       nddmspc   I    number of DDM species
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c       10/24/01  Removed BSWAP and converted integer strings to character*4
c       07/19/02  Added seperate source area map for each grids.
c       05/01/03  Time span of emissions must now match emiss update interval
c       07/16/07 -bkoo-     Added HRVOC
c       03/01/16 -gwilson-  Added partial source area map
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
c       
c-----------------------------------------------------------------------
c   Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      integer ndate
      real    ttime
      integer nox
      integer noy
      integer nlay_ems
      integer nddmspc
      real    emisddm(nox,noy,nlay_ems,nddmspc)
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  cname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, igroup, idx, iseg
      integer       i, j, nedge, iptr, ispc, iddm, izone, nlays_in
      integer       imap, imod, ipart, idxfile, buffer_offset
      integer       nx, ny
      real          emstmp, bgtim, edtim, orgx, orgy, utmx, utmy
      real          dx, dy
      logical       lfound, luse, lpass
c
      real emsgrd(MXCELLS,MXCELLS)
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c   --- initialize emissions to zero ---
c
      emisddm = 0.
c
c  --- if not doing emissions groups, just return here ---
c
      if( nemddm .EQ. 0 ) goto 9999
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
      else
        nedge = 1
      endif
c
c   --- loop over all of the groups ----
c
      do 10 igroup=1,ngroup
       do 11 idxfile=1,num_iortem(igrid,igroup)
c
c   --- skip if filename not supplied ---
c
        if( is_netcdf_iortem(igrid,igroup,idxfile) .OR. 
     &      .NOT. ltemfl(igrid,igroup,idxfile) .OR. .NOT. larsrc ) cycle
c
c   --- set the unit number for surface emissions file ---
c
        iounit = iortem(igrid,igroup,idxfile)
        fname = temfil(igrid,igroup,idxfile)
        buffer_offset = buffer_offset_iortem(igrid,igroup,idxfile)
c
c   --- get the number of layers for this file ---
c
        rewind(iounit)
        read(iounit)
        read(iounit) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nlays_in
        read(iounit)
        read(iounit)
c
c   --- read the date and time, again ---
c
        lfound = .FALSE.
        lpass = .FALSE.
  111   continue
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
            write(iout,'(//,a)')'ERROR in READARDDM:'
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
     &            ((emsgrd(i,j),i=1,nox-2*buffer_offset),j=1,noy-2*buffer_offset)
             write(cname,'(10A1)') iname
             do k=1,nlays_in
                do i=1,nox
                   emsgrd(i,1) = emsgrd(i,1)
                   emsgrd(i,noy) = emsgrd(i,noy-2*buffer_offset)
                enddo
                do j=1,noy
                   emsgrd(1,j) = emsgrd(1,j)
                   emsgrd(nox,j) = emsgrd(nox-2*buffer_offset,j)
                enddo
             enddo
c
c   --- if date and time does not match this hour, skip this record ---
c
             if( le1day ) then
                 if( bgtim .NE. ttime ) goto 20
             else
                 if( ndate .NE. ibgdat .OR. bgtim .NE. ttime ) goto 20
             endif
             lfound  = .TRUE.
c
c   --- if the species is a not modelled or not used for DDM ----
c
             idx = idxems(igrid,igroup,idxfile,ispc)
             if( idx .LE. 0 ) goto 20
c
c  --- find this species in the modeled species list ---
c
             luse = .FALSE.
             do imod=1,nspec
                if( cname .EQ. spname(imod) ) then
                   do iddm=1,nemddm
                      luse = .FALSE.
                      if( emddmsp(iddm) .EQ. cname ) luse = .TRUE.
                      if( emddmsp(iddm) .EQ. NAMALL  ) luse = .TRUE.
                      if( emddmsp(iddm) .EQ. NAMVOC
     &                               .AND. lvocsp(imod) ) luse = .TRUE.
                      if( emddmsp(iddm) .EQ. NAMNOX
     &                               .AND. lnoxsp(imod) ) luse = .TRUE.
                     if( emddmsp(iddm) .EQ. NAMHRV
     &                               .AND. lhrvoc(imod) ) luse = .TRUE.
c
c  --- if this DDM species matches this modeled species, load it ---
c
                      if( luse ) then
                         do 30 j=2,noy-1
                            do 40 i=2,nox-1
                               imap = igrmap(0,1,igrid,i,j)
                               if( imap .LE. 0 .OR. imap .GT. nregin )
     &                                                          goto 40
c
c  --- convert to emissions time and put into array ---
c
                               emstmp = emsgrd(i,j)/(60.*dtems)
                               iptr = iptddm(imod) + nicddm +
     &                             nbcddm*nedge + (iddm-1)*ngroup*nregin +
     &                                      (imap-1)*ngroup + igroup - 1
                               emisddm(i,j,1,iptr) = emisddm(i,j,1,iptr) + emstmp
c
c  --- next cell ---
c
   40                       continue
   30                    continue
                      endif
c
c  --- next DDM emssions species ---
c
                   enddo
c
c  --- next modeled species ---
c
                endif
             enddo
c
c   --- next species ---
c
  20        continue
c
c   --- if the correct hour has not been found, 
c       go back and read some more else read next file ---
c
            if( .NOT. lfound ) then
               goto 111
            else
               do i=1,nspcem(igrid,igroup,idxfile) + 1
                   backspace(iounit)
               enddo
               goto 11
            endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222       continue
            if( le1day ) then
                rewind(iounit)
                read(iounit)
                read(iounit)
                read(iounit)
                read(iounit)
                if( lpass ) goto 7001
                if( .NOT. lpass ) lpass = .TRUE.
                goto 111
               else
                goto 7001
            endif
c
c   --- next file in this group/grid ---
c
  11     continue
c
c  --- get the next group ---
c
  10   continue
       goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in RDARDDM:' 
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'ERROR: Reading emissions ',
     &                  'after hour ',ibgdat, bgtim,' in file: ',fname
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in RDARDDM:' 
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 
     &              'ERROR: Premature end-of-file',
     &              ' in emissions file after hour ',ibgdat, bgtim,
     &                                               ' in file: ',fname
      call camxerr()
c
c----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
