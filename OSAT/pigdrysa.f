c*** PIGDRYSA
c
      subroutine pigdrysa(ncolx,nrowy,nlays,nrad,nspc,ntrac,
     &                    icell,jcell,idxpig,
     &                    deldry,saconc,dryfld)
      use grid
      use pigsty
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine puts the dry deposition from the PiG into the tracer
c     dry deposition array.
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
c       idxpig  I  the index of the puff in PiG arrays
c       deldry  R  change in model dry depostion
c       saconc  R  gridded array of tracer concentrations
c       dryfld  R  gridded array of tracer dry depositions
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c       02/10/16  --bkoo--      Created
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
      integer   idxpig
      real      deldry(nspc)
      real      saconc(ncolx,nrowy,nlays,ntrac)
      real      dryfld(ncolx,nrowy,ntrac)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer icls, ispc, itrc, i
      integer iOO, iO3
      real    delcls(MXALCLS)
      real    sumcls
      real    sumO3, sumOO, prodOO, destOO
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- get dry dep change for model species of each tracer class
c
      do icls = 1, ntrcls
         delcls(icls) = 0.
         do ispc = nrad+1, nspc
            delcls(icls) = delcls(icls)
     &                   + deldry(ispc) * fluxmap(ispc,icls)
         enddo
      enddo
c
c  --- update dry dep for O3 & OO tracers
c
      if ( lozone ) then

         sumO3 = 0.
         do i = iptcls(idxipt(ITRO3N)), nptcls(idxipt(ITRO3V))
            sumO3 = sumO3 + saconc(icell,jcell,1,i)
         enddo
         sumO3 = MAX(BNDLPT, sumO3)

         sumOO = 0.
         do i = iptcls(idxipt(ITROON)), nptcls(idxipt(ITROOV))
            sumOO = sumOO + saconc(icell,jcell,1,i)
         enddo
         sumOO = MAX(BNDLPT, sumOO)

         prodOO = MAX( 0.,  delcls(idxipt(ITRRGN)) )
         destOO = MAX( 0., -delcls(idxipt(ITRRGN)) )

         if ( delcls(idxipt(ITRO3V)) .GT. 0. ) then
c
c  --- this is not supposed to happen, but just in case ...
c      If this happens, just allocate it to O3V of this plume
c
            itrc = iemcls(idxipt(ITRO3V)) - 1 + 
     &                    ipufmap(idxpig) + ipufgrp(idxpig) * nregin

            dryfld(icell,jcell,itrc) = dryfld(icell,jcell,itrc)
     &                               + delcls(idxipt(ITRO3V))
         else
c
c  --- allocate O3 decrease based on relative contribution of O3X
c
            do i = iptcls(idxipt(ITRO3N)), nptcls(idxipt(ITRO3V))
               dryfld(icell,jcell,i) = dryfld(icell,jcell,i)
     &                               + delcls(idxipt(ITRO3V))
     &                               * saconc(icell,jcell,1,i) / sumO3
            enddo
         endif
c
c  --- allocate OO increase based on relative contribution of O3X
c           and OO decrease based on relative contribution of OOX
c
         iO3 = iptcls(idxipt(ITRO3N)) - 1
         do i = iptcls(idxipt(ITROON)), nptcls(idxipt(ITROOV))
            iO3 = iO3 + 1
            dryfld(icell,jcell,i) = dryfld(icell,jcell,i)
     &                    + prodOO * saconc(icell,jcell,1,iO3) / sumO3
     &                      - destOO * saconc(icell,jcell,1,i) / sumOO
         enddo
c
c  --- reset the dry dep change for the tracers already processed so that they will be skipped below
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
         if ( delcls(icls) .EQ. 0. ) CYCLE
c
c  --- for all other tracers
c
         if ( delcls(icls) .GT. 0. ) then
c
c  --- if positive change, allocate it to this plume
c
            itrc = iemcls(icls) - 1 + 
     &                    ipufmap(idxpig) + ipufgrp(idxpig) * nregin

            dryfld(icell,jcell,itrc) = dryfld(icell,jcell,itrc)
     &                               + delcls(icls)
         else
c
c  --- if negative change, reduce dry dep based on relative contribution of the tracer
c
            sumcls = 0.
            do i = iptcls(icls), nptcls(icls)
               sumcls = sumcls + saconc(icell,jcell,1,i)
            enddo
            sumcls = MAX(BNDLPT, sumcls)

            do i = iptcls(icls), nptcls(icls)
               dryfld(icell,jcell,i) = dryfld(icell,jcell,i)
     &                              + delcls(icls)
     &                              * saconc(icell,jcell,1,i) / sumcls
            enddo
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
