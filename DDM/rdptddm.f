c***** RDPTDDM.F
c
      subroutine rdptddm(ndate,ttime)
      use filunit
      use grid
      use chmstry
      use ptemiss
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
c     Inputs:
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     10/24/01  Removed BSWAP and converted integer strings to character*4
c     07/19/02  Added seperate source area map for each grids.
c     05/01/03  Time span of emissions must now match emiss update interval
c     07/16/07 -bkoo-     Added HRVOC
c     03/01/16 -gwilson-  Added partial source area map
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
      integer       ibgdat, iendat, iounit, igroup, idx, iseg
      integer       i, j, nedge, iptr, ispc, iddm, igrd
      integer       icel, jcel, imap, npoint, idum, imod, izcel(MXPTSRC)
      integer       idcompact(MXPTSRC), idxfile
      real          bgtim, edtim, emstmp, pi
      real          xloccrs, yloccrs, xloctmp, yloctmp
      logical       lfound, luse, lpass
c
      real emspnt(MXPTSRC), effph_in(MXPTSRC), flowrat_in(MXPTSRC)
c
      common /ddmptdat/ emspnt, effph_in, flowrat_in
c
      data pi /3.1415927/
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c  --- initialize emissions to zero ---
c
      do i=1,nptsrc
         do j=1,ntotsp
            sapnts(i,j) = 0.
         enddo
      enddo
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
        do 11 idxfile=1,num_iortpt(igroup)
c
c   --- skip if filename not supplied ---
c
          if( is_netcdf_iortpt(igroup,idxfile) .OR. .NOT. 
     &            ltptfl(igroup,idxfile) .OR. .NOT. lptsrc ) cycle
c
c   --- set the unit number for surface emissions file ---
c
          iounit = iortpt(igroup,idxfile)
          fname = tptfil(igroup,idxfile)
c
c   --- read the date and time, again ---
c
          lfound = .FALSE.
          lpass = .FALSE.
  111     continue
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
              write(iout,'(//,a)')'ERROR in RDPTDDM:'
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
          read(iounit,ERR=7000,END=7001) iseg, npoint
c
          if( npoint .GT. MXPTSRC ) goto 7002
          read(iounit,ERR=7000,END=7001) (idcompact(i), idum,
     &                izcel(i), flowrat_in(i), effph_in(i),i=1,npoint)
          do i=1,npoint
            idxcp = i
            if( lcompactpt(igroup,idxfile) ) idxcp = idcompact(i)
            if(flowrat_in(idxcp) .LT. 0. .AND. effph_in(idxcp) .LT. 0. ) then
              if( ABS(effph_in(idxcp)) .LT. ABS(flowrat_in(idxcp)) ) then
                 write(iout,'(//,a)')'ERROR in READPT:'
                 write(iout,*) 'Invalid values found for flow rate and plume height.'
                 write(iout,*) 'The plume distribution override is triggered. But the '
                 write(iout,*) 'value for flow rate (base) is larger than value for ',
     &                      'plume height (top).'
                 call camxerr()
              endif
            endif
            idxbase = idx_point_in_list(igroup,idxfile,idxcp)
            effph(idxbase) = effph_in(idxcp)
            flowrat(idxbase) = flowrat_in(idxcp)
          enddo
c 
          do 20 ispc=1,nspcpt(igroup,idxfile)
              read(iounit,ERR=7000) iseg, (iname(i),i=1,10), 
     &                                        (emspnt(i),i=1,npoint)
              write(cname,'(10A1)') iname
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
               idx = idxpts(igroup,idxfile,ispc)
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
     &                                 .AND. lvocsp(imod) ) luse = .TRUE.
                        if( emddmsp(iddm) .EQ. NAMNOX
     &                                 .AND. lnoxsp(imod) ) luse = .TRUE.
                        if( emddmsp(iddm) .EQ. NAMHRV
     &                                 .AND. lhrvoc(imod) ) luse = .TRUE.
c
c  --- if this DDM species matches this modeled species, load it ---
c
                        if( luse ) then
                           do 30 i=1,npoint
                           if( emspnt(i) .LE. 0 ) cycle
                              idxcp = i
                              if( lcompactpt(igroup,idxfile) ) idxcp = idcompact(i)
                              idxbase = idx_point_in_list(igroup,idxfile,idxcp)
                              if( llatlon ) then
                                 xloccrs = xlocpt(idxbase) - xorg
                                 yloccrs = ylocpt(idxbase) - yorg
                              else
                                 xloccrs = xlocpt(idxbase)/1000. - xorg
                                 yloccrs = ylocpt(idxbase)/1000. - yorg
                              endif
                              icel = 1 + FLOOR( xloccrs/delx )
                              jcel = 1 + FLOOR( yloccrs/dely )
                              if( icel .LT. 1 .OR. icel .GT. ncol(1) ) 
     &                                                            goto 30
                              if( jcel .LT. 1 .OR. jcel .GT. nrow(1) ) 
     &                                                          goto 30
c
c   --- find out if a nest contains this source  ---
c
                              igrd = 1
                              do ig = 2,ngrid
                                 xloctmp = xloccrs - (inst1(ig)-1)*delx
                                 yloctmp = yloccrs - (jnst1(ig)-1)*dely
                                 if( xloctmp .GT. 0. .AND. 
     &                                           yloctmp .GT. 0. ) then
                                      ii = 2 + FLOOR( xloctmp/delx *
     &                                              FLOAT( meshold(ig) ) )
                                      jj = 2 + FLOOR( yloctmp/dely *
     &                                              FLOAT( meshold(ig) ) )
                                      if( ii .GT. 1 .AND. jj .GT. 1 .AND. 
     &                                     ii .LT. ncol(ig) .AND. 
     &                                             jj .LT. nrow(ig) ) then
                                         igrd = ig
                                         icel = ii
                                         jcel = jj
                                      endif
                                 endif
                              enddo
c
c  --- get the region for this cell from mapping array ----
c
                              imap = igrmap(0,1,igrd,icel,jcel)
                              if( izcel(i) .LT. 0 ) then
                                 if( .NOT. lptoverride ) goto 7003
                                 imap = ABS( izcel(i) )
                              endif
                              if( imap .LE. 0 .OR. imap .GT. nregin ) 
     &                                                         goto 30
c
c  --- convert to emissions time and put into array ---
c
                              emstmp = emspnt(i)/(60.*dtems)
                              iptr = iptddm(imod) + nicddm + 
     &                            nbcddm*nedge + (iddm-1)*ngroup*nregin +
     &                                        (imap-1)*ngroup + igroup - 1
                                if( emstmp .GT. bdnl(imod) )
     &                                 sapnts(idxbase,iptr) = 
     &                                     sapnts(idxbase,iptr) + emstmp
c
c  --- next point ---
c
   30                      continue
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
  20       continue
c
c   --- if the correct hour has not been found, 
c       go back and read some more else read next file ---
c
           if( .NOT. lfound ) then
              goto 111
           else
              do i=1,nspcpt(igroup,idxfile) + 3
                  backspace(iounit)
              enddo
              goto 11
           endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222      continue
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
c  --- get the next file ---
c
  11     continue
  10  continue
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)')'ERROR in RDPTDDM:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'ERROR: Reading emissions ',
     &    'after hour ',ibgdat, bgtim,' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)')'ERROR in RDPTDDM:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)')
     &              'ERROR: Premature end-of-file',
     &              ' in point source file after hour ',ibgdat, bgtim,
     &                                ' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDPTDDM:'
      write(iout,'(/,1X,A,I10,2A)') 'Number of points: ',npoint,
     &                   ' exceeds max in file: ',fname(:istrln(fname))
      write(iout,'(1X,2A,/,1X,A)') 'Increase parameter MXPTSRC.'
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in RDPTDDM:'
      write(iout,'(/,1X,2A )') 'A source was read that has the point ',
     &                   'source override trigger turned on.'
      write(iout,'(1X,2A)') 'If you intend to use point source ',
     &                     'override please set the namelist variable'
      write(iout,'(1X,2A)') 'DDM_PT_Override to TRUE in the CAMx_Control namelist.'
      write(iout,'(1X,2A)') 'If you do not want point source override ',
     &                         'you need to set the kcell variable '
      write(iout,'(1X,2A)') 'in your point source file to zero.'
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
