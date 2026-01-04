c*** PIGDUMPSA
c
      subroutine pigdumpsa(ncolx,nrowy,nlays,nrad,nspc,ntrac,
     &                     icell,jcell,kclbeg,kclend,idxpig,
     &                     delconc,saconc)
      use grid
      use pigsty
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine puts the mass from the PiG into the tracer
c     concentration array.  The mass is added to the tracer species
c     from the region/group from which the PiG originated.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c   Argument descriptions:
c     Inputs:
c       ncolx   I  number of columns
c       nrowy   I  number of rows
c       nlays   I  number of layers
c       nrad    I  number of radical species
c       nspc    I  number of model species
c       ntrac   I  number of tracers
c       icell   I  the X grid location of current cell
c       jcell   I  the X grid location of current cell
c       kclbeg  I  the bottom vertical layer in the affected column
c       kclend  I  the top vertical layer in the affected column
c       idxpig  I  the index of the puff in PiG arrays
c       delconc R  change in model concentrations 
c       saconc  R  gridded array of tracer concentrations
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c       12/08/96  --gwilson--   Original development
c       11/10/97  --gwilson--   Removed unused argument: IGRID
c       08/25/05  --cemery--    Revamped to use puff-specific source 
c                               region/group pointers
c       02/10/16  --bkoo--      Updated for OSAT3
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   ncolx
      integer   nrowy
      integer   nlays
      integer   nrad
      integer   nspc
      integer   ntrac
      integer   icell
      integer   jcell
      integer   kclbeg
      integer   kclend
      integer   idxpig
      real      delconc(MXLAYER,MXSPEC+1)
      real      saconc(ncolx,nrowy,nlays,ntrac)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer icls, ispc, itrc, i, kcel
      integer iOO, iO3
      real    delcls(MXALCLS)
      real    sumcls
      real    sumO3, sumOO, prodOO, destOO, OOtoO3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- loop over all layers in the affected column ---
c  
      do kcel = kclbeg, kclend
c
c  --- get concentration change for model species of each tracer class
c
         do icls = 1, ntrcls
            delcls(icls) = 0.
            do ispc = nrad+1, nspc
               delcls(icls) = delcls(icls)
     &                      + delconc(kcel,ispc) * fluxmap(ispc,icls)
            enddo
         enddo
c
c  --- update O3 & OO tracers
c
         if ( lozone ) then

            sumO3 = 0.
            do i = iptcls(idxipt(ITRO3N)), nptcls(idxipt(ITRO3V))
               sumO3 = sumO3 + saconc(icell,jcell,kcel,i)
            enddo
            sumO3 = MAX(BNDLPT, sumO3)

            sumOO = 0.
            do i = iptcls(idxipt(ITROON)), nptcls(idxipt(ITROOV))
               sumOO = sumOO + saconc(icell,jcell,kcel,i)
            enddo
            sumOO = MAX(BNDLPT, sumOO)

            prodOO = MAX( 0.,  delcls(idxipt(ITRRGN)) )
            destOO = MAX( 0., -delcls(idxipt(ITRRGN)) )

            if ( delcls(idxipt(ITRO3V)) .GT. 0. ) then
c
c  --- this is not supposed to happen, but just in case ...
c
c  --- allocate O3 from NO2 destruction based on relative contribution of OOX
c
               OOtoO3 = MIN( destOO, delcls(idxipt(ITRO3V)) )

               iOO = iptcls(idxipt(ITROON)) - 1
               do i = iptcls(idxipt(ITRO3N)), nptcls(idxipt(ITRO3V))
                  iOO = iOO + 1
                  saconc(icell,jcell,kcel,i) = saconc(icell,jcell,kcel,i)
     &                  + OOtoO3 * saconc(icell,jcell,kcel,iOO) / sumOO
               enddo
c
c  --- allocate any remaining O3 increase to O3V of this plume
c
               itrc = iemcls(idxipt(ITRO3V)) - 1 + 
     &                       ipufmap(idxpig) + ipufgrp(idxpig) * nregin

               saconc(icell,jcell,kcel,itrc) = saconc(icell,jcell,kcel,itrc)
     &                                + delcls(idxipt(ITRO3V)) - OOtoO3
            else
c
c  --- allocate O3 decrease based on relative contribution of O3X
c
               do i = iptcls(idxipt(ITRO3N)), nptcls(idxipt(ITRO3V))
                  saconc(icell,jcell,kcel,i) = saconc(icell,jcell,kcel,i)
     &                             + delcls(idxipt(ITRO3V))
     &                             * saconc(icell,jcell,kcel,i) / sumO3
                  saconc(icell,jcell,kcel,i) =
     &                          MAX(BNDLPT, saconc(icell,jcell,kcel,i))
               enddo
            endif
c
c  --- allocate OO increase based on relative contribution of O3X
c           and OO decrease based on relative contribution of OOX
c
            iO3 = iptcls(idxipt(ITRO3N)) - 1
            do i = iptcls(idxipt(ITROON)), nptcls(idxipt(ITROOV))
               iO3 = iO3 + 1
               saconc(icell,jcell,kcel,i) = saconc(icell,jcell,kcel,i)
     &                  + prodOO * saconc(icell,jcell,kcel,iO3) / sumO3
     &                    - destOO * saconc(icell,jcell,kcel,i) / sumOO
               saconc(icell,jcell,kcel,i) =
     &                          MAX(BNDLPT, saconc(icell,jcell,kcel,i))
            enddo
c
c  --- reset the conc change for the tracers already processed so that they will be skipped below
c
            delcls(idxipt(ITRO3N)) = 0.
            delcls(idxipt(ITRO3V)) = 0.
            delcls(idxipt(ITROON)) = 0.
            delcls(idxipt(ITROOV)) = 0.
         endif
c
c  --- loop over all tracer classes
c
         do icls = 1, ntrcls
c
c  --- skip if no change or already processed
c
            if ( ABS(delcls(icls)) .LT. BNDLPT ) CYCLE
c
c  --- for all other tracers
c
            if ( delcls(icls) .GT. 0. ) then
c
c  --- if positive change, allocate it to this plume
c
               itrc = iemcls(icls) - 1 + 
     &                       ipufmap(idxpig) + ipufgrp(idxpig) * nregin

               saconc(icell,jcell,kcel,itrc) = saconc(icell,jcell,kcel,itrc)
     &                                       + delcls(icls)
            else
c
c  --- if negative change, reduce the tracer based on its relative contribution
c
               sumcls = 0.
               do i = iptcls(icls), nptcls(icls)
                  sumcls = sumcls + saconc(icell,jcell,kcel,i)
               enddo
               sumcls = MAX(BNDLPT, sumcls)

               do i = iptcls(icls), nptcls(icls)
                  saconc(icell,jcell,kcel,i) = saconc(icell,jcell,kcel,i)
     &                            + delcls(icls)
     &                            * saconc(icell,jcell,kcel,i) / sumcls
                  saconc(icell,jcell,kcel,i) =
     &                          MAX(BNDLPT, saconc(icell,jcell,kcel,i))
               enddo
            endif
c
c  --- next tracer class ---
c
         enddo
c
c  --- next affected layer ---
c
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
