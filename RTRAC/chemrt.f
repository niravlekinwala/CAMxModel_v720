c*** CHEMRT
c
      subroutine chemrt(numcol,numrow,numlay,igrid,icell,jcell,kcell, 
     &                  i0,j0,pcell,tcell,cold,cnew,oh,o3,no3,dtime,
     &                  irt_cel)
      use chmstry
      use camxcom
      use grid
      use rtracchm
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine performs chemistry adjustments on the reactive
c     tracer gasses using information from the host chemistry solution.  
c
c      Copyright 1996 - 2022
c     Ramboll
c
c   Argument descriptions:
c     Inputs:
c       numcol  I number of columns in this slice
c       numrow  I number of rows in this slice
c       numlay  I number of layers in this slice
c       igrid   I  number of the grid containing the cell
c       icell   I  X index of the cell
c       jcell   I  Y index of the cell
c       kcell   I  Z index of the cell
c       pcell   R  cell pressure (mb)
c       tcell   R  cell temperature (K)
c       cold    R  species concentrations before host chemistry (umol/m3)
c       cnew    R  species concentrations after host chemistry (umol/m3)
c       oh      R  OH radical concentration (ppm)
c       o3      R  ozone concentration (ppm)
c       no3     R  no3 radical concentration (ppm)
c       dtime   R  time step (hours)
c       irt_cel I  index for RTRAC receptor cells
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     01/16/02   --gwilson--    Original development
c     11/12/03   --cemery--     Added secondary RTRAC species to primary
c                               RTRAC species
c     09/22/08   --gyarwood--   Delete CONVFAC from calculation of RTRAC
c                               species secondary to a regular species
c     11/06/12   --gwilson--    Fixed receptors for MPI
c     04/05/21   --cemery--     Input reaction rates are molecular units
c                               and converted to ppm units here
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
      integer numcol
      integer numrow
      integer numlay
      integer igrid
      integer icell
      integer jcell
      integer kcell
      real    pcell
      real    tcell
      real    cold(*)
      real    cnew(*)
      real    oh
      real    o3
      real    no3
      real    dtime
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   i
      integer*8 idx
      integer*8 idxcel
      real      factor
      real      prod
      real      rkoh
      real      rko3
      real      rkno3
      real      rkspc
c
      real drate(MXTRSP)
      real rtnew(MXTRSP)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- calculate the index of the cell in the grid ---
c
      idxcel = DBLE(ipsa3d(igrid)-1) + DBLE(icell) + DBLE(numcol) * 
     &          DBLE(jcell-1) + DBLE(numcol) * DBLE(numrow) * DBLE(kcell-1)
      do i=1,nrtrac
         idx = idxcel + DBLE(numcol) * DBLE(numrow) *
     &                                      DBLE(numlay) * DBLE(i-1)
         ptconc(idx) = max(ptconc(idx),rtlbnd(i))
      enddo
c
c --- factor is a temp/press adjustment for second order
c     rate constants that convert from cm3 molec-1 sec-1 to
c     ppm-1 min-1, also convert to per hour ---
c
      factor = 60.*pcell*100./(8.314*tcell)*6.022E+23/1.0E12 * 60.
c
c --- calculate the total decay rate and any secondary production
c     rate for each tracer.  get host photolysis rates from rk in the
c     chmstry.inc include file ---
c
      do i=1,nrtgas
         if( jnum(i) .GT. 0 ) then
            rkspc = rk(jnum(i))
         else
            rkspc = 0.0
         endif 
         rkoh  = arrhen(aoh(i),eaoh(i),boh(i),troh(i),tcell) * factor
         rko3  = arrhen(ao3(i),eao3(i),bo3(i),tro3(i),tcell) * factor
         rkno3 = arrhen(ano3(i),eano3(i),bno3(i),trno3(i),tcell)
     &           * factor
         drate(i) = oh*rkoh + o3*rko3 + no3*rkno3 + rkspc*rtjfact(i)
         idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                      DBLE(numlay) * DBLE(i-1)
         rtnew(i) = ptconc(idx)*EXP(-drate(i)*dtime)
         rtnew(i) = MAX(rtlbnd(i),rtnew(i))
         prod = 0.0
         if( lsecnd(i) )then
            if( lreg(i) ) then
               prod = cnew(ksec(i)) - cold(ksec(i))*EXP(-drate(i)*dtime)
               prod = MAX(prod,0.0)
            else
               idk = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(ksec(i)-1)
               prod = AMAX1(ptconc(idk) - rtnew(ksec(i)),0.0)
            endif
         endif
         rtnew(i) = rtnew(i) + prod
      enddo
c
      do i = 1,nrtgas
         idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
         ptconc(idx) = rtnew(i)
      enddo
c
c --- if this is a receptor cell, accumulate average decay rates
c
      if( lrcpfil .AND. irt_cel .GT. 0 ) then 
         do i=1,nrtgas
            rcpdcy(irt_cel, i) = rcpdcy(irt_cel,i) + 
     &                              drate(i)*dtime/(dtout/60.)
         enddo
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
c
      function arrhen(a,ea,b,tref,temp)
c
c --- Calculate a rate constant using k = A*(T/Tref)^B*exp(-Ea/T)
c
      real arrhen
c
      arrhen = a
      if( tref .NE. 0. ) arrhen = a*((temp/tref)**b)*exp(-ea/temp)
c
      end
