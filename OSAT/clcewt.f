c*** CLCEWT
c
      subroutine clcewt(jdate,etim)
      use filunit
      use grid
      use chmstry
      use bndary
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine calculates the weighted reactivity factor for VOC
c   species in each of the source groups.  All of the emissions for the
c   group are read and the emissions are weighted by reactivity factor
c   and summed up.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Argument description:
c        jdate   I   ending date of simulation (YYJJJ)
c        etim    R   ending time of simulation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/04/96   --gwilson--    Original development
c     10/10/96   --gwilson--    Added code to output emissions and 
c                               reactivity for each grouping.
c     12/08/96   --gwilson--    Added code to set the index into tracer
c                               species list for the PiG sources
c     12/12/96   --gwilson--    Fixed bug in reporting total NOx Tons
c     01/09/97   --gwilson--    Fixed (another) bug in reporting total 
c                               NOx Tons
c     01/12/97   --gwilson--    Fixed bug in calculating the source
c                               region in the fine grid
c     11/06/01   --cemery--     Input dates are now Julian
c     10/28/09   --gwilson--    Made the large local arrays
c                               allocatable to avoid memory issues
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   jdate
      real      etim
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c  CVTTON   R   conversion factor for grams to tons
c
      real   CVTTON
c
      parameter( CVTTON = 907184.7 )
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer idx
      integer i, j, ndlast, ncola, nrowa
      integer ncount, ioff, icls, itrc
      logical lerror
      real    sumtrc, sumkoh, summir, yldhvoc, yldlvoc, sumcls, difmax, diff
      real    ttlast, sumyld, allmax
c
      real    tonems(MXTRSP), tolerance_val

      real, allocatable, dimension(:,:) :: emssum
      real, allocatable, dimension(:,:) :: emsbas
      real, allocatable, dimension(:,:) :: emsoth

      real*8, allocatable, dimension(:,:,:) :: emslft
      real*8, allocatable, dimension(:,:,:) :: emstot

      logical lemit(MXALCLS)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      tolerance_val = 0.0017
c
c  --- allocate the local arrays ---
c
      allocate( emssum(nspec,nsaspc) )
      allocate( emsbas(nspec,nsaspc) )
      allocate( emsoth(nspec,nsaspc) )
c
      ncola = maxval( ncol(1:ngrid) )
      nrowa = maxval( nrow(1:ngrid) )
      allocate( emslft(ncola,nrowa,nspec) )
      allocate( emstot(ncola,nrowa,nspec) )
c
c   --- set the date and times ---
c
      lerror = .FALSE.
      ndlast = jdate 
      ttlast = etim/100.0
      if( ttlast .EQ. 0.0 ) then
         ttlast = 24.0
         ndlast = ndlast - 1 
      endif
c
c   --- initialize the array to zero ---
c
      tonems = 0.
      emssum = 0.
      emsbas = 0.
      emsoth = 0.
      emstot = 0.
      emslft = 0.
      lemit = .FALSE.
c
c   --- loop over all of the groups ----
c
      call sumgrps(ncola,nrowa,nspec,nsaspc,ndlast,ttlast*100.,emstot,
     &                           emslft,emsbas,emsoth,emssum,lemit)
      call ncf_sumgrps(ncola,nrowa,nspec,nsaspc,ndlast,ttlast,emstot,
     &                           emslft,emsbas,emsoth,emssum,lemit)
c
c  --- check that all emissions are accounted for (could have some
c      machine fuzz) ----
c
      allmax = -99999999.
      do idx=1,nspec
         difmax = -99999999.
         do j=2,nrow(1)-1
            do i=2,ncol(1)-1
               if( emstot(i,j,idx) .NE. 0. ) then 
                  diff = ABS( (emstot(i,j,idx) - 
     &                        emslft(i,j,idx)) / emstot(i,j,idx) )
                  if( diff .GT. difmax ) then
                      imax = i
                      jmax = j
                      difmax = diff
                  endif
                  allmax = MAX(diff,allmax)
               endif
            enddo
         enddo
         if( .NOT. leftovr_area .AND. (difmax .GT. tolerance_val)
     &                                       .AND. .NOT. lerror ) then
             write(iout,'(//,a)') 'ERROR in CLCEWT:'
             if( emslft(imax,jmax,idx) .LT. emstot(imax,jmax,idx) ) then
                 write(iout,'(/,1X,2A)') 'There is a significant amount ',
     &             'of emissions unaccounted for in source groupings.'
                 write(iout,'(1X,2A)') 'You should turn on the ',
     &            'Use_Gridded_Leftover_Group group flag in job script.'
             else
                 write(iout,'(/,1X,2A)') 'There is an inconsistency between ',
     &             'model emissions and source group emissions.'
                 write(iout,'(1X,2A)') 'Please check all of your emissions files.'
             endif
             write(iout,'(1X,2A)') 'See the .diag file for an ',
     &                                             'emissions table.'
             write(iout,'(/,A,A)') 'Species: ',spname(idx)
             write(iout,'(A,2I5)') 'Cell Index: ',imax,jmax
             write(iout,'(A,E11.5)') 'Value in regular model files  : ',
     &                                            emstot(imax,jmax,idx) 
             write(iout,'(A,E11.5)') 'Value in emissions group files: ',
     &                                            emslft(imax,jmax,idx) 
             lerror = .TRUE.
         endif
      enddo
      if( leftovr_area .AND. (allmax .LE. 0.0005) .AND. .NOT. lerror ) then
          write(iout,'(//,a)') 'ERROR in CLCEWT:'
          write(iout,'(/,1X,2A)') 'The "leftover" emissions group ',
     &                       'has an insignificant amount of emissions.'
          write(iout,'(1X,2A)') 'You should turn off the ',
     &                        '"leftover_area group" flag in job script.'
          write(iout,'(1X,2A)') 'See the .diag file for an ',
     &                                                'emissions table.'
          lerror = .TRUE.
      endif
c
c  --- all emissions are summed, calculate the weghted fraction ----
c
      do icls=1,ntrcls
        do 10 i=iptcls(icls),nptcls(icls)
c
c   --- ignore if this is an initial condition or boundary condition
c       tracer ---
c
           if( ptname(i)(7:8) .EQ. 'IC' ) goto 10
           if( ptname(i)(7:8) .EQ. 'BC' ) goto 10
c
c   --- find the class for this species ---
c
           sumtrc = 0.
           sumkoh = 0.
           summir = 0.
           sumyld = 0.
           yldhvoc = 0.
           yldlvoc = 0.
           do idx=1,nspec
              if( emssum(idx,i) .GT. 0 ) then
                 sumtrc = sumtrc + emssum(idx,i) * trspmap(idx,icls)
                 sumkoh = sumkoh + emssum(idx,i) * rkohrt(idx) * 
     &                                           trspmap(idx,icls)
                 summir = summir + emssum(idx,i) * rmirrt(idx) * 
     &                                           trspmap(idx,icls)
                 if( yhratmap(idx,icls) .GT. 0. .OR.
     &               ylratmap(idx,icls) .GT. 0. ) then
                    sumyld = sumyld + emssum(idx,i)
                    yldhvoc = yldhvoc + emssum(idx,i) * yhratmap(idx,icls)
                    yldlvoc = yldlvoc + emssum(idx,i) * ylratmap(idx,icls)
                 endif
                 if( trspmap(idx,icls) .GT. 0. ) tonems(i) = 
     &                  tonems(i) + emssum(idx,i) * mwspec(idx) / CVTTON
              endif
           enddo
           if( sumtrc .GT. 0. ) then
              wtkoh(i) = sumkoh / sumtrc
              wtmir(i) = summir / sumtrc
           else
              wtkoh(i) = 0.
              wtmir(i) = 0.
           endif
           if( sumyld .GT. 0. ) then
              yhrates(i) = yldhvoc / sumyld
              ylrates(i) = yldlvoc / sumyld
           else
              yhrates(i) = 0.
              ylrates(i) = 0.
           endif
  10    continue
      enddo
c
c  --- calculate the "leftover" group from lump sums ---
c
      if( leftovr_area ) then
         do icls=1,ntrcls
            do i=1,nregin
               sumtrc = 0.
               sumkoh = 0.
               summir = 0.
               sumyld = 0.
               yldhvoc = 0.
               yldlvoc = 0.
               do idx=1,nspec
                  itrc = iemcls(icls) - 1 + i + ngroup*nregin
                  diff = emsbas(idx,itrc) - emsoth(idx,itrc)
                  if( diff .GT. 0. ) then
                      sumtrc = sumtrc + diff * trspmap(idx,icls)
                      sumkoh = sumkoh + diff * rkohrt(idx) * 
     &                                            trspmap(idx,icls)
                      summir = summir + diff * rmirrt(idx) * 
     &                                            trspmap(idx,icls)
                      if( yhratmap(idx,icls) .GT. 0. .OR.
     &                    ylratmap(idx,icls) .GT. 0. ) then
                         sumyld = sumyld + diff
                         yldhvoc = yldhvoc + diff * yhratmap(idx,icls)
                         yldlvoc = yldlvoc + diff * ylratmap(idx,icls)
                      endif
                      if( trspmap(idx,icls) .GT. 0. ) tonems(itrc) = 
     &                      tonems(itrc) + diff * mwspec(idx) / CVTTON
                  endif
                enddo
                if( sumtrc .GT. 0 ) then
                    wtkoh(itrc) = sumkoh / sumtrc
                    wtmir(itrc) = summir / sumtrc
                else
                    wtkoh(itrc) = 0.
                    wtmir(itrc) = 0.
                endif
                if( sumyld .GT. 0 ) then
                    yhrates(itrc) = yldhvoc / sumyld
                    ylrates(itrc) = yldlvoc / sumyld
                else
                    yhrates(itrc) = 0.
                    ylrates(itrc) = 0.
                endif
            enddo
         enddo
      endif
c
c  --- echo the data if doing ozone ---
c
      do icls=1,ntrcls
          if( idxcls(icls) .EQ. ITROON .OR.
     &        idxcls(icls) .EQ. ITROOV ) CYCLE 
          if( ngroup .EQ. 0 ) then
              ioff = 0
              ncount = 0
          else
              if( leftovr_area ) then
                  ncount = ngroup + 1
              else
                  ncount = ngroup
              endif
              ioff = 1
          endif
          do i=ioff,ncount
             if( lemit(icls) ) then
                write(idiag,9000) 'Species   ','Group ','Region ',
     &                  '     Average ','Emissions' 
                write(idiag,9001)' ','Reactivity','(Tons)'
                write(idiag,9002) ('-',j=1,60)
                sumcls = 0.
                do j=1,nregin
                   if( i .GT. 0 ) then
                      itrc = iemcls(icls) - 1 + j + (i-1)*nregin
                   else
                      itrc = iemcls(icls) - 1 + j
                   endif
                   if( tonems(itrc) .LE. 0. ) then
                      write(idiag,9003) ptname(itrc),i,j,
     &                                       wtkoh(itrc),tonems(itrc)
                   else if( tonems(itrc) .GT. 0. .AND.
     &                                   tonems(itrc) .LE. 0.01 ) then
                      write(idiag,9006) ptname(itrc),i,j,
     &                                       wtkoh(itrc),tonems(itrc)
                   else if( tonems(itrc) .GT. 0.01 .AND.
     &                                     tonems(itrc) .LE. 10. ) then
                      write(idiag,9004) ptname(itrc),i,j,
     &                                       wtkoh(itrc),tonems(itrc)
                   else if( tonems(itrc) .GT. 10. .AND. 
     &                                  tonems(itrc) .LE. 999999.) then
                      write(idiag,9005) ptname(itrc),i,j,
     &                                       wtkoh(itrc),tonems(itrc)
                   else 
                      write(idiag,9006) ptname(itrc),i,j,
     &                                       wtkoh(itrc),tonems(itrc)
                   endif
                   sumcls = sumcls + tonems(itrc)
                enddo
                write(idiag,9002) ('-',j=1,60)
                if( sumcls .LE. 0. ) then
                   write(idiag,9013) 'Total',0.,sumcls
                else if( sumcls .GT. 0. .AND. sumcls .LE. 10. ) then
                   write(idiag,9014) 'Total',0.,sumcls
                else if( sumcls .GT. 10. .AND. sumcls .LE. 999999.) then
                   write(idiag,9015) 'Total',0.,sumcls
                else
                   write(idiag,9016) 'Total',0.,sumcls
                endif
             endif
          enddo
       enddo
       write(idiag,'(//)')
c
c  --- deallocate the local arrays ---
c
      deallocate( emssum )
      deallocate( emsbas )
      deallocate( emsoth )
c
      deallocate( emslft )
      deallocate( emstot )
c
c  --- return to the calling routine ---
c
      if( lerror ) call camxerr()
      goto 9999
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,1X,A10,3X,A,3X,A,3X,A13,4X,A)
 9001 format(1X,A10,25X,A10,3X,A10,3X,A10,3X,A10)
 9002 format(100(A1))
 9003 format(1X,A10,4X,I3,5X,I3,6X,F10.1,3X,F10.0)
 9004 format(1X,A10,4X,I3,5X,I3,6X,F10.1,3X,F10.3)
 9005 format(1X,A10,4X,I3,5X,I3,6X,F10.1,3X,F10.1)
 9006 format(1X,A10,4X,I3,5X,I3,6X,F10.1,3X,E10.4)
 9013 format(1X,A10,21X,F10.1,3X,F10.0)
 9014 format(1X,A10,21X,F10.1,3X,F10.3)
 9015 format(1X,A10,21X,F10.1,3X,F10.1)
 9016 format(1X,A10,21X,F10.1,3X,E10.4)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
