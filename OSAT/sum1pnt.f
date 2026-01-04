      subroutine sum1pnt(numcols,numrows,nspmod,nsptrac,igroup,idxfile,
     &             idx,numpts,emsbas,emsoth,emslft,emstot,
     &                             emspnt,emssum,izcel,idcompact,lemit)
      use grid
      use chmstry
      use pigsty
      use filunit
      use tracer
c
c      Copyright 1996 - 2022
c     Ramboll
c
c
c----CAMx v7.20 220430
c
c     SUM1PNT sums up the point emission of one species for a given group
c
c       07/19/02  --gwilson-- Added seperate source area map for each grid.
c       08/25/05  --cemery--  Revamped PiG pointer arrays for source group
c                             and region
c       11/16/06  --gwilson-- fixed bug in point source override for PiG sources
c       11/27/06  --gwilson-- fixed bug in calculating emissions table
c       10/28/09  --gwilson-- Changed dimension of variables to accomodate 
c                             the dynamic memory allocation
c       12/20/13  --gwilson-- Added compact point source file
c       03/01/16  --gwilson-- Added partial source area map
c       01/17/22  --gwilson-- Added conditional to bypass overrided 
c                             check when using only 1 region
c
c     Input argument:
c        numcols           max number of columns in any grid
c        numrows           max number of columns in any grid
c        nspmod            number of model species
c        nsptrac           number of tracer species
c        igroup            group ID
c        idxfile           file number
c        idx               specie ID
c        emspnt            the species emission for the group
c        idcompact         index of the source in master list
c
c     Output arguments:
c        emssum            emission summed over grid
c        emslft            leftover emission
c        emstot            total emission
c        lemit             flag to determine if tracer class is emitted
c
      include "camx.prm"
      include "flags.inc"
c
      character*200 fname
      integer   numcols
      integer   numrows
      integer   nspmod
      integer   nsptrac
      integer   igroup
      integer   idxfile
      integer   idx
      integer   numpts
      real*8    emslft(numcols,numrows,nspmod)
      real*8    emstot(numcols,numrows,nspmod)
      real      emspnt(MXPTSRC)
      real      emsbas(nspmod,nsptrac)
      real      emsoth(nspmod,nsptrac)
      real      emssum(nspmod,nsptrac)
      integer   izcel(*)
      integer   idcompact(*)
      logical   lemit(*)
      logical   luse
      integer   ipart
      real      frac
      integer   istart
c
      real, allocatable, dimension(:) :: xloctmp
      real, allocatable, dimension(:) :: yloctmp
c
c  --- allocate the local arrays ---
c
      allocate( xloctmp(numpts) )
      allocate( yloctmp(numpts) )
c
c   --- set the file name for elevated points emissions file ---
c
      num_ptsfiles = num_iortpt(igroup)
      fname = tptfil(igroup,idxfile)
c
c  --- make sure this species is needed ---
c
      luse = .FALSE.
      do icls=1,ntrcls
        if( trspmap(idx,icls) .NE. 0.  .OR.
     &                         yhratmap(idx,icls) .NE. 0. .OR.
     &                        ylratmap(idx,icls) .NE. 0. ) then
          if( trspmap(idx,icls) .NE. 0. ) lemit(icls) = .TRUE.
          luse = .TRUE.
        endif
      enddo
      if( .NOT. luse ) goto 111
c
c   --- sum up the emissions for each point ---
c
      istart = idx_start_sapts(igroup,idxfile)
      if (llatlon) then
        do n = istart+1,istart+nptsrc_safile(igroup,idxfile)
          idxpt = idcompact(n)
          xloctmp(n) = xlocpt(n) - xorg
          yloctmp(n) = ylocpt(n) - yorg
        enddo
      else
        do n = istart+1,istart+nptsrc_safile(igroup,idxfile)
          idxpt = idcompact(n)
          xloctmp(n) = xlocpt(n)/1000. - xorg
          yloctmp(n) = ylocpt(n)/1000. - yorg
        enddo
      endif
      n = 0
      do 70 i=istart+1,istart+nptsrc_safile(igroup,idxfile)
         n = n+1
         idxpt = idcompact(n)
         icel = 1 + FLOOR( xloctmp(i)/delx )
         jcel = 1 + FLOOR( yloctmp(i)/dely )
         if(icel .LE. 0 .OR. icel .GT. ncol(1)) goto 70
         if(jcel .LE. 0 .OR. jcel .GT. nrow(1)) goto 70
         icrs = icel
         jcrs = jcel
c
c   --- find out if a nest contains this source  ---
c
         igrd = 1
         do ig = 2,ngrid
           xlocnst = xloctmp(i) - (inst1(ig)-1)*delx
           ylocnst = yloctmp(i) - (jnst1(ig)-1)*dely
           ii = 2 + FLOOR( xlocnst/delx * FLOAT( meshold(ig) ) )
           jj = 2 + FLOOR( ylocnst/dely * FLOAT( meshold(ig) ) )
           if( ii .GT. 1 .AND. jj .GT. 1 .AND. ii .LT. ncol(ig) .AND.
     &                                           jj .LT. nrow(ig) ) then
              igrd = ig
              icel = ii
              jcel = jj
            endif
         enddo
c
c  --- get the region for this cell from mapping array ----
c
        imap = igrmap(0,1,igrd,icel,jcel)
        frac = 1.0
c
c  --- change the region if the override is set ---
c
        if( nregin .GT. 1 ) then
           if( izcel(idcompact(n)) .LT. 0 ) then
               if( .NOT. lptoverride ) goto 7000
               imap = ABS( izcel(idcompact(n)) )
               frac = 1.0
           endif
        endif
        if( imap .LE. 0 .OR. imap .GT. nregin ) goto 7001
c
c  --- calculate the index into the tracer species for this group/region ---
c
        if( ngroup .GT. 0 ) then
c
c  --- if doing PiG and source is a PIG source, set the PiG map and group
c      pointers ---
c
          if( lpigsa(i) .AND. emspnt(n) .GT. 0. ) then
             ipigmap(i) = imap
             ipiggrp(i) = igroup-1
          endif
c
          do icls=1,ntrcls
             if( trspmap(idx,icls) .NE. 0.  .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                if( leftovr_area ) then
                    ipt = iemcls(icls)-1 + imap + ngroup*nregin
                    emsbas(idx,ipt) = emsbas(idx,ipt) + emspnt(n) * frac
                    emsoth(idx,ipt) = emsoth(idx,ipt) + emspnt(n) * frac
                endif
                ipt = iemcls(icls)-1 + imap +(igroup-1)*nregin
                emssum(idx,ipt) = emssum(idx,ipt) + emspnt(n) * frac
             endif
          enddo
          do icls=1,ntrcls
              if( trspmap(idx,icls) .NE. 0.  .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                emslft(icrs,jcrs,idx) = 
     &                         emslft(icrs,jcrs,idx) + emspnt(n) * frac
                emstot(icrs,jcrs,idx) = 
     &                         emstot(icrs,jcrs,idx) + emspnt(n) * frac
             endif
          enddo
c
c   --- only using regular model emissions ---
c
        else
           do icls=1,ntrcls
              if( trspmap(idx,icls) .NE. 0.  .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                 ipt = iemcls(icls) - 1 + imap
                 emssum(idx,ipt) = emssum(idx,ipt) + emspnt(n) * frac
              endif
           enddo
c
c  --- if doing PiG and source is a PIG source, set the PiG map and group
c      pointers ---
c
           if( lpigsa(n) .AND. emspnt(n) .GT. 0. ) then
              ipigmap(n) = imap
              ipiggrp(n) = 0
           endif
        endif
  70  continue
c
c  --- processed entire file ---
c
 111  continue
c
c  --- deallocate the local arrays ---
c
      deallocate( xloctmp )
      deallocate( yloctmp )
      goto 999

c
c----------------------------------------------------------------------
c  Error messages:
c----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in SUM1PNT:'
      write(iout,'(/,1X,2A )') 'A source was read that has the point ',
     &                   'source override trigger turned on.'
      write(iout,'(1X,2A)') 'If you intend to use point source ',
     &                     'override please set the namelist variable'
      write(iout,'(1X,2A)') 'SA_PT_Override to TRUE in the CAMx_Control namelist.'
      write(iout,'(1X,2A)') 'If you do not want point source override ',
     &                         'you need to set the kcell variable '
      write(iout,'(1X,2A)') 'in your point source file to zero.'
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in SUM1PNT:'
      write(iout,'(/,1X,A,A,I4)') 'Invalid region found in',
     &             ' point source override when reading point source file.'
      write(iout,'(1X,A,I4)') 'Region code      : ',imap
      write(iout,'(1X,A,I4)') 'Number of regions: ',nregin
      write(iout,'(10X,A,/,A)') 'Point source filename: ',
     &                                                 fname(:istrln(fname))
      write(iout,'(1X,2A)') 'Check the values in the point ',
     &                                                     'source overide.'
      call camxerr()
c
 999  continue
      return
      end
