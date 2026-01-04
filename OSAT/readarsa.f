c*** READARSA.F
c
      subroutine readarsa(igrid,idate,btim,nox,noy,nlay_ems,
     &                                            nspsa,dx,dy,emisar)
      use filunit
      use grid
      use bndary
      use camxcom
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions from each of the surface
c   emissions files and calculates the NOx and VOC levels.  It then
c   places these emissions in the appropriate place in the gridded
c   array used for tracer emissions.  The tracer to which the NOx and
c   VOC emissions is assigned depends on the source group and the
c   source region.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Argument declarations:
c        Outputs:
c          emisar    R  array to store the tracer emissions
c        Inputs:
c          igrid     I  grid number for this grid
c          idate     I  date of current hour (YYJJJ)
c          btim      R  time of current hour
c          nox       I  number of X cells in grid
c          noy       I  number of Y cells in grid
c          nlays_ems I  number of layers in emissions files
c          nspsa     I  number of Y cells in grid
c          dx        R  cell width in X direction
c          dy        R  cell width in Y direction
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     06/06/96   --gwilson--    Original development
c     01/09/97   --gwilson--    Fixed a bug in fineding the source
c                               region in the fine grid
c     01/12/97   --gwilson--    Now checks for negative emissions in
c                               leftover group and exits if found
c     02/01/97   --gwilson--    Put fuzz factor of 0.1 ppb to determine
c                               if emissions are truly negative.
c     02/03/97   --gwilson--    Put code to ignore emissions in the
c     11/06/01   --cemery--     Input dates are now Julian
c     07/19/02   --gwilson--    Seperate source area map for each grids
c     12/04/02   --gyarwood--   Improved message for negative leftover
c     03/01/16   --gwilson--    Added partial source area map
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
      integer   igrid
      integer   idate
      real      btim
      integer   nox
      integer   noy
      integer   nlay_ems
      integer   nspsa
      real*4    dx(noy)
      real*4    dy
      real      emisar(nox,noy,nlay_ems,nspsa)
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   FUZZ  R  fuzz factor used to determine if emissions are truly < 0.
c
      real   FUZZ
c
      parameter( FUZZ = -0.0001 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      ndate, imap, itim, ipart
      integer      i, j, k, l
      real         ttime, frac
c
      real, allocatable, dimension(:,:,:,:) :: lefcls
      real, allocatable, dimension(:,:,:,:) :: emscls
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- allocate the temporary space ---
c
      allocate( lefcls(ntrcls,nox,noy,nlay_ems) )
      allocate( emscls(ntrcls,nox,noy,nlay_ems) )
c
c   --- set the date and times ---
c
      ndate = idate
      ttime = btim/100.0
c
c   --- if up to time for releasing a new timing tracer, bump up
c       counter for number of tracer species ---
c
      if( igrid .EQ. 1 .AND. ntrtim .GT. 0 ) then
          if( MOD( INT(ttime), INT(24/ntrtim) ) .EQ. 0 ) then
              nreles = nreles + 1
              nsaspc = nsaspc + 2 * nregin
          endif
      endif
c
c   --- initialize emissions to zero ---
c
      emisar = 0.
c
c  --- for each group, call routine to read the emissions for
c      the group ----
c
      do igroup = 0,ngroup
         call rdargrp(igrid,ndate,ttime*100.,nox,noy,nlay_ems,
     &                                       ntrcls,igroup,emscls)
         call ncf_rdargrp(igrid,ndate,ttime*100,nox,noy,nlay_ems,
     &                                       ntrcls,igroup,emscls)
c
c  ---  if doing source groups calculate the left-over and put emissions 
c       into proper place in the emissions  array ---
c
         if( ngroup .GT. 0 ) then
            do 10 j=2,noy-1
               do 20 i=2,nox-1
c
c  --- get the region for this cell from mapping array ----
c
                  do ipart=1,npartial(igroup,igrid)
                     imap = igrmap(igroup,ipart,igrid,i,j)
                     frac = frcmap(igroup,ipart,igrid,i,j)
                     if( (imap .LE. 0 .AND. frac .GT. 0.) .OR. imap .GT. nregin ) goto 20
c
c  --- calculate the leftover emissions to use in last source group ---
c
                     if( igroup .EQ. 0 ) then
                        do k=1,nlay_ems
                           do icls=1,ntrcls
                              lefcls(icls,i,j,k) = emscls(icls,i,j,k) * frac
                           enddo
                        enddo
                     else
                        do k=1,nlay_ems
                           do icls=1,ntrcls
                              lefcls(icls,i,j,k) = lefcls(icls,i,j,k) -
     &                                              emscls(icls,i,j,k) * frac
c
c  --- put emissions into position in gridded tracer emissions array ---
c
                              ipos = iemcls(icls)-1 + imap+(igroup-1)*nregin
                              emisar(i,j,k,ipos) = emisar(i,j,k,ipos) + 
     &                                             emscls(icls,i,j,k) * frac
                           enddo
                        enddo
                     endif
                  enddo
   20          continue
   10        continue
c
c  --- only one group, just load emissions into arrays from
c      postion 0 in gridded array ----
c
         else
             do 30 j=2,noy-1
                 do 40 i=2,nox-1
c
c  --- get the region for this cell from mapping array ----
c
                    do ipart=1,npartial(igroup,igrid)
                        imap = igrmap(0,ipart,igrid,i,j)
                        frac = frcmap(0,ipart,igrid,i,j)
                        if( (imap .LE. 0 .AND. frac .GT. 0.) .OR. imap .GT. nregin ) goto 40
c
c   --- put emissions in array at correct offset ---
c
                        do k=1,nlay_ems
                           do icls=1,ntrcls
                              ipos = iemcls(icls) - 1 + imap
                              emisar(i,j,k,ipos) = emisar(i,j,k,ipos ) + 
     &                                             emscls(icls,i,j,k) * frac
                            enddo
                        enddo
                     enddo
   40            continue
   30        continue
          endif
c
c  --- read emissions for the next group ---
c
      enddo
c
c  --- put leftover emissions in last group ----
c
      if( ngroup .GT. 0 .AND. leftovr_area ) then
          do 50 j=2,noy-1
            do 60 i=2,nox-1
c
c  --- get the region for this cell from mapping array ----
c
                imap = igrmap(0,1,igrid,i,j)
                frac = frcmap(0,1,igrid,i,j)
                if( (imap .LE. 0 .AND. frac .GT. 0.) .OR. imap .GT. nregin ) goto 60
c
c   ---- make sure leftover group is not negative ---
c
                do k=1,nlay_ems
                   do icls=1,ntrcls
                      if( lefcls(icls,i,j,k) .LT. 0. ) then
                        if( lefcls(icls,i,j,k) .GT. FUZZ ) then
                            lefcls(icls,i,j,k) = 0.
                         else
                            write(iout,'(//,A)') 'ERROR in READARSA:'
                            write(iout,'(/,4A,I3,A,I3,A,I3,2A,I3)')
     &                             'Negative ',clsnam(idxcls(icls)),
     &                          ' emissions calculated in leftover group ',
     &                                        'in cell: (',i,',',j,',',k,') ',
     &                                         'in Grid: ', igrid
                            write(iout,'(A,F20.6)')'   Value = ',
     &                                           lefcls(icls,i,j,k)
                            write(iout,'(A,I2,A,F20.6)') ' group ',
     &                                  igroup,' = ', emscls(icls,i,j,k)
                            call camxerr()
                         endif
                      endif
                      ipos = iemcls(icls) - 1 + imap + ngroup * nregin
                      emisar(i,j,k,ipos) = emisar(i,j,k,ipos ) + 
     &                                       lefcls(icls,i,j,k) * frac
                   enddo
                enddo
   60        continue
   50     continue
       endif
c
c   --- deallocate the temporary space ---
c
      deallocate( lefcls )
      deallocate( emscls )
c
      return
      end
