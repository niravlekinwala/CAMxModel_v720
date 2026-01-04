c*** PIGWETSA
c
      subroutine pigwetsa(ncolx,nrowy,nlays,nrad,nspc,ntrac,
     &                    icell,jcell,kclbeg,kclend,idxpig,
     &                    delwet,saconc,wetfld)
      use grid
      use pigsty
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine puts the wet deposition from the PiG into the tracer
c     wet deposition array.
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
c       delwet  R  change in model wet depostion
c       saconc  R  gridded array of tracer concentrations
c       wetfld  R  gridded array of tracer wet depositions
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
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
      real      delwet(nspc)
      real      saconc(ncolx,nrowy,nlays,ntrac)
      real      wetfld(ncolx,nrowy,ntrac)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer icls, ispc, itrc, i, kcel
      integer iO3
      real    delcls(MXALCLS)
      real    sumO3, fracO3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- get wet dep change for model species of each tracer class
c
c      Note: IRONWET enforces DELWET >= 0, therefore always DELCLS >= 0
c
      do icls = 1, ntrcls
         delcls(icls) = 0.
         do ispc = nrad+1, nspc
            delcls(icls) = delcls(icls)
     &                   + delwet(ispc) * fluxmap(ispc,icls)
         enddo
      enddo
c
c  --- update wet dep for O3 & OO tracers
c
      if ( lozone ) then
c
c  --- O3 is not supposed to increase, but just in case ...
c      If this happens, just allocate it to O3V of this plume
c
         itrc = iemcls(idxipt(ITRO3V)) - 1 + 
     &                 ipufmap(idxpig) + ipufgrp(idxpig) * nregin

         wetfld(icell,jcell,itrc) = wetfld(icell,jcell,itrc)
     &                            + delcls(idxipt(ITRO3V))
c
c  --- allocate OO increase based on relative contribution of O3X
c      across the layers this puff covers (layer weighting is ignored though)
c      assuming that RGN increase is completely due to NOx titration of O3
c      (this sounds inconsistent because no wet dep is allowed for negative O3 increment)
c
         sumO3 = 0.
         do i = iptcls(idxipt(ITRO3N)), nptcls(idxipt(ITRO3V))
            do kcel = kclbeg, kclend
               sumO3 = sumO3 + saconc(icell,jcell,kcel,i)
            enddo
         enddo
         sumO3 = MAX(BNDLPT, sumO3)

         iO3 = iptcls(idxipt(ITRO3N)) - 1
         do i = iptcls(idxipt(ITROON)), nptcls(idxipt(ITROOV))
            iO3 = iO3 + 1
            fracO3 = 0.
            do kcel = kclbeg, kclend
               fracO3 = fracO3 + saconc(icell,jcell,kcel,iO3) / sumO3
            enddo
            wetfld(icell,jcell,i) = wetfld(icell,jcell,i)
     &                            + delcls(idxipt(ITRRGN)) * fracO3
         enddo
c
c  --- reset the wet dep change for the tracers already processed so that they will be skipped below
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
         if ( delcls(icls) .GT. 0. ) then
c
c  --- allocate positive change to this plume
c
            itrc = iemcls(icls) - 1 + 
     &                    ipufmap(idxpig) + ipufgrp(idxpig) * nregin

            wetfld(icell,jcell,itrc) = wetfld(icell,jcell,itrc)
     &                               + delcls(icls)
         endif
c
c  --- next tracer class ---
c
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
