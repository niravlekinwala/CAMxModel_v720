      subroutine sum1grd(ncola,nrowa,numcols,numrows,numlays,
     &                  nspmod,nsptrac,igroup,igrd,idx,emssum,
     &                    emsgrd,emsbas,emsoth,emslft,emstot,lemit)
      use grid
      use chmstry
      use tracer
c
c
c      Copyright 1996 - 2022
c     Ramboll
c
c
c----CAMx v7.20 220430
c
c     SUM1GRD sums up the area emission of one species for a given group
c     in a given grid
c
c       07/19/02  --gwilson-- Added seperate source area map for each grids.
c       10/28/09  --gwilson-- Changed dimension of variables to accomodate
c                             the dynamic memory allocation
c       03/01/16  --gwilson-- Added partial source area map
c
c     Input argument:
c        numcols           number of columns in any grid
c        numrows           number of columns in any grid
c        numlays           number of columns in any grid
c        nspmod            number of model species
c        nsptrac           number of tracer species
c        igroup            group ID
c        igrd              grid ID
c        idx               specie ID
c        emsgrd            the species emission in the grid
c
c     Output arguments:
c        emssum            emission summed over grid
c        emsbas            base emission
c        emsoth            "otherwise" emission
c        emslft            leftover emission
c        emstot            total emission
c        lemit             flag to determine if tracer class is emitted
c        
      include "camx.prm"
c
      integer   ncola
      integer   nrowa
      integer   numcols
      integer   numrows
      integer   numlays
      integer   nspmod
      integer   nsptrac
      real      emssum(nspmod,nsptrac)
      real      emsgrd(numcols,numrows,numlays)
      real      emsbas(nspmod,nsptrac)
      real      emsoth(nspmod,nsptrac)
      real*8    emslft(ncola,nrowa,nspmod)
      real*8    emstot(ncola,nrowa,nspmod)
      logical   lemit(*)
      integer   ipart
      real      frac
c
      logical   luse
c
c  --- Entry point ---
c
c  --- make sure this species is used ---
c
      luse = .FALSE.
      do icls=1,ntrcls
        if( trspmap(idx,icls) .NE. 0.  .OR.
     &                         yhratmap(idx,icls) .NE. 0. .OR.
     &                         ylratmap(idx,icls) .NE. 0. ) then
          if( trspmap(idx,icls) .NE. 0. ) lemit(icls) = .TRUE. 
          luse = .TRUE.
        endif
      enddo
      if( .NOT. luse ) goto 9999
c
c   --- sum up the emissions excluding the boundary for this grid ---
c
      do 40 j=2,nrow(igrd)-1
        jcrs = (j - 2)/nmesh(igrd) + j1(igrd)
        if( igrd .eq. 1 ) jcrs = j
        do 50 i=2,ncol(igrd)-1
c
c   --- calculate the coarse grid offset ---
c
           icrs = (i - 2)/nmesh(igrd) + i1(igrd)
           if( igrd .eq. 1 ) icrs = i
c
c   --- skip cell if this grid has a child in this cell ---
c
           ijcl = i + (j-1)*ncol(igrd)
           if( idfin(iptr2d(igrd)-1+ijcl) .NE. 0 ) goto 50
c
c  --- get the region for this cell from mapping array,
c      the grid cell should be the coarse grid  ----
c
           do ipart=1,npartial(igroup,igrd)
              imap = igrmap(igroup,ipart,igrd,i,j)
              frac = frcmap(igroup,ipart,igrd,i,j)
              if( (imap .LE. 0 .AND. frac .GT. 0.) .OR. imap .GT. nregin ) goto 50
c
c  --- calculate the index into the tracer species for this gruoup/region ---
c
              if( ngroup .GT. 0 ) then
c
c   --- if group is base emissions, add to "leftover" group ----
c
                 if( igroup .EQ. 0 ) then
                   if( leftovr_area ) then
                      do icls=1,ntrcls
                        if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                           ipt = iemcls(icls)-1 + imap+ngroup*nregin
                           do k=1,numlays
                             emsbas(idx,ipt) = emsbas(idx,ipt) + emsgrd(i,j,k) * frac
                           enddo
                        endif
                      enddo
                   endif
                   do icls=1,ntrcls
                     if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                       do k=1,numlays
                         emstot(icrs,jcrs,idx) = 
     &                              emstot(icrs,jcrs,idx) + emsgrd(i,j,k) * frac
                       enddo
                     endif
                   enddo
c
c   --- otherwise, add to this group/region and subtract from "leftover" ---
c
                 else
                   do icls=1,ntrcls
                     if( trspmap(idx,icls) .NE. 0. .OR. 
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                       ipt = iemcls(icls)-1 + imap+(igroup-1)*nregin
                       do k=1,numlays
                           emssum(idx,ipt) = emssum(idx,ipt) + emsgrd(i,j,k) * frac
                           if( leftovr_area ) then
                              ipt = iemcls(icls)-1 + imap+ngroup*nregin
                              emsoth(idx,ipt) = emsoth(idx,ipt) + emsgrd(i,j,k) * frac
                           endif
                        enddo
                     endif
                   enddo
                   do icls=1,ntrcls
                     if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                       do k=1,numlays
                          emslft(icrs,jcrs,idx) = 
     &                               emslft(icrs,jcrs,idx) + emsgrd(i,j,k) * frac
                       enddo
                     endif
                enddo
              endif
c
c   --- only using regular model emissions ---
c
           else
                 do icls=1,ntrcls
                   if( trspmap(idx,icls) .NE. 0. .OR.
     &                                 yhratmap(idx,icls) .NE. 0. .OR.
     &                                 ylratmap(idx,icls) .NE. 0. ) then
                      ipt = iemcls(icls) - 1 + imap
                      do k=1,numlays
                        emssum(idx,ipt) = emssum(idx,ipt) + emsgrd(i,j,k) * frac
                      enddo
                   endif
                 enddo
              endif
           enddo
  50    continue
  40  continue
c
 9999 continue
      return
      end
