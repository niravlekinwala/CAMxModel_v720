c***** RDPTRT.F
c
      subroutine rdptrt(ndate,ttime)
      use filunit
      use grid
      use chmstry
      use ptemiss
      use bndary
      use camxcom
      use rtracchm
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
c   This routine reads one hour of emissions for the RTRAC process
c   and fills the approproate arrays.  The emissions file for one grid
c   but each emissions groups is read.
c    Argument descriptions:
c     Outputs:
c     Inputs:
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/24/01   --gwilson--   Originial development
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
      integer ndate
      real    ttime
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  cname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, idx, iseg
      integer       i, j, ispc
      integer       icel, jcel, npoint, idum, imod
      real          bgtim, edtim, xloctmp, yloctmp
      real          emstmp
      logical       lfound, lpass
c
      real emspnt(MXPTSRC), effph_in(MXPTSRC), flowrat_in(MXPTSRC)
c
      common /rtptdat/ emspnt, effph_in, flowrat_in
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c   --- skip if filename not supplied ---
c
      if( is_netcdf_iortpt(1,1) .OR. .NOT. ltptfl(1,1) 
     &                          .OR. .NOT. lptsrc ) goto 9999
c
c  --- initialize emissions to zero ---
c
      do i = 1,nptsrc
        do j = 1,ntotsp
         sapnts(i,j) = 0.
        enddo
      enddo
c
c   --- set the unit number for surface emissions file ---
c
      iounit = iortpt(1,1)
      fname = tptfil(1,1)
c
c   --- read the date and time, again ---
c
      lfound = .FALSE.
      lpass = .FALSE.
  111 continue
      read(iounit,END=222) ibgdat, bgtim, iendat, edtim
      bgtim = 100.*aint(bgtim) + 60.*amod(bgtim,1.)
      edtim = 100.*aint(edtim) + 60.*amod(edtim,1.)
c
c   --- read the emissions for this hour ---
c

      read(iounit,ERR=7000,END=7001) iseg, npoint
      if( npoint .GT. nptsrc ) goto 7002
      read(iounit,ERR=7000,END=7001) (idum, idum,
     &                idum, flowrat_in(i), effph_in(i),i=1,npoint)
      do i=1,npoint
        if(flowrat_in(i) .LT. 0. .AND. effph_in(i) .LT. 0. ) then
           if( ABS(effph_in(i)) .LT. ABS(flowrat_in(i)) ) then
              write(iout,'(//,a)')'ERROR in READPT:'
              write(iout,*) 'Invalid values found for flow rate and plume height.'
              write(iout,*) 'The plume distribution override is triggered. But the '
              write(iout,*) 'value for flow rate (base) is larger than value for ',
     &                   'plume height (top).'
              call camxerr()
           endif
        endif
        idxbase = idx_point_in_list(1,1,i)
        effph(idxbase) = effph_in(i)
        flowrat(idxbase) = flowrat_in(i)
      enddo
c
      do 20 ispc = 1,nspcpt(1,1)
        read(iounit,ERR=7000) iseg, (iname(i),i=1,10), 
     &                              (emspnt(i),i=1,npoint)
        write(cname,'(10A1)') iname
c
c   --- if date and time does not match this hour, skip this record ---
c
        if( le1day ) then
           if( abs(bgtim-ttime).ge.0.01 ) goto 20
        else
           if( ndate .NE. ibgdat .OR. abs(bgtim-ttime).ge.0.01 ) goto 20
        endif
        lfound  = .TRUE.
c
c   --- if the species is not modelled or not used for RTRAC ----
c
        idx = idxpts(1,1,ispc)
        if( idx .LE. 0 ) goto 20
c
c  --- find this species in the modeled species list ---
c
        do imod = 1,nrtrac
          if( cname .EQ. ptname(imod) ) then
             do 30 i = 1,npoint
               idxbase = idx_point_in_list(1,1,i)
               if( llatlon ) then
                  xloctmp = xlocpt(i) - xorg
                  yloctmp = ylocpt(i) - yorg
               else
                  xloctmp = xlocpt(i)/1000. - xorg
                  yloctmp = ylocpt(i)/1000. - yorg
               endif
               icel = 1 + FLOOR( xloctmp/delx )
               jcel = 1 + FLOOR( yloctmp/dely )
               if( icel .LT. 1 .OR. icel .GT. ncol(1) ) goto 30
               if( jcel .LT. 1 .OR. jcel .GT. nrow(1) ) goto 30
c
c  --- convert to emissions time and put into array ---
c
               emstmp = emspnt(i)/(60.*dtems)
               if( emstmp .GT. rtlbnd(imod) ) sapnts(idxbase,imod) = emstmp
c
c  --- next point ---
c
   30        continue
             goto 20
          endif
c
c  --- next RTRAC emssions species ---
c
        enddo
c
c  --- if not found in model species list, write a message ----     
c   
        write(idiag,'(1X,4A)') 'Species in RTRAC ',
     &            'point source file: ',cname(:istrln(cname)),
     &            ' not found in species list ... Skipping.'
c
c   --- next species ---
c
  20  continue
c
c   --- if the correct hour has not been found, 
c       go back and read some more else read next file ---
c
      if( .NOT. lfound ) then
         goto 111
      else
         goto 9999
      endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222 continue
      if( le1day ) then
          rewind(iounit)
          read(iounit)
          read(iounit)
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
c  --- get the next file ---
c
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in RDPTRT:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'ERROR: Reading emissions ',
     &    'after hour ',ibgdat, bgtim,' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in RDPTRT:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)')
     &              'ERROR: Premature end-of-file',
     &              ' in point source file after hour ',ibgdat, bgtim,
     &                                ' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDPTGRP:'
      write(iout,'(/,1X,A,I10,2A)') 'Number of points: ',npoint,
     &     ' is not consistent with regular emissions in file: ',
     &                                            fname(:istrln(fname))
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
