      subroutine partitionrt(numcol,numrow,numlay,nspc,igrd,icell,jcell,
     &                       kcell,conc)
      use filunit
      use chmstry
      use grid
      use rtracchm
      use tracer
c
c----CAMx v7.20 220430
c 
c-----------------------------------------------------------------------
c   Description:
c     This routine performs gas-aerosol partitioning of RTRAC gas species 
c     using Koa theory.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c   Argument descriptions:
c     Inputs:
c       numcol  I number of columns in this slice
c       numrow  I number of rows in this slice
c       numlay  I number of layers in this slice
c       nspc    I number of host model species
c       igrd    I  umber of the grid containing the cell
c       icell   I  X index of the cell
c       jcell   I  Y index of the cell
c       kcell   I  Z index of the cell
c       conc    R  host model species concentrations (umol/m3, ug/m3)
c
c     Routines called: 
c             
c     Called by: 
c        CHEMDRIV
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     10/14/09  --jjung--   Original development
c      9/18/15  --cemery--  Bug fix on PTCONC gas units
c      8/25/16  --bkoo--    Updated for new SOAP
c      5/12/21  --cemery--  Minor update to skip species with Koa = 0
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c 
      implicit none
c     include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer numcol
      integer numrow
      integer numlay
      integer nspc
      integer igrd
      integer icell
      integer jcell
      integer kcell
      real    conc(nspc)
      real    eps
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   i
      integer   ispc
      integer*8 idx
      integer*8 idxaer
      integer*8 idxcel

      real daero
      real Mom,fom,Kp,TSP
      real MwH2O,Mw,gas,aero
      real total
c
      data MwH2O /18.0/   !Mw of H2O (g/mol)
      data eps /1.0e-20/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- Calculate array pointers for the current grid cell ---
c
      idxcel =  DBLE(ipsa3d(igrd)-1) + 
     &                 DBLE(icell) + DBLE(numcol)*DBLE(jcell-1) + 
     &                       DBLE(numcol)*DBLE(numrow)*DBLE(kcell-1)
c
c --- Calculate organic and total PM mass, and OM/PM ratio
c
      TSP = 0.0
      do ispc = ngas+1,nspc
         TSP = TSP + conc(ispc)
      enddo
      TSP = max(TSP,eps)
      if (TSP.eq.eps) then
         fom = 0.2 ! Set 20%
      else
         Mom = conc(ksoa1) + conc(ksoa2) + conc(ksoa3) + conc(ksoa4) +
     &         conc(ksopa) + conc(ksopb) + conc(kpoa)
         fom = Mom/TSP
      endif
      fom = max(fom,1.e-5)        ! Bound fom to 0.001%
c
c-----Loop over GAS species, and calculate partitioning for this cell
c
      do i = 1,nrtgas
         idx    = idxcel + DBLE(numcol)*DBLE(numrow)*DBLE(numlay)*DBLE(i-1)
         idxaer = idxcel + DBLE(numcol)*DBLE(numrow)*DBLE(numlay)*DBLE(nrtgas+i-1)
         ptconc(idx) = max(ptconc(idx),rtlbnd(i))
         ptconc(idxaer) = max(ptconc(idxaer),rtlbnd(nrtgas+i-i))
         if (eqkoa(i).eq.0.) cycle
         Mw = rtmolwt(i)
         gas = ptconc(idx)*Mw                   ! umol/m3 -> ug/m3
         total = gas + ptconc(idxaer)
         Kp = 1.0e-9/820.0*eqkoa(i)*fom
         aero = total*Kp*TSP/(Kp*TSP + 1.)
         daero = aero - ptconc(idxaer)
         if (daero.ge.0.0) then                         ! aerosol increase
            daero = min(daero,gas-rtlbnd(i)*Mw)
            if (daero.ge.gas) then                      ! Negative gas
              write(iout,*)'Error in PARTITIONRT:',
     &                     ' Partitioning more than available gas.'
              write(iout,*)'Available gas (ug/m3)',gas
              write(iout,*)'Mass subtracted (ug/m3)',daero
              write(iout,*)'Kp,Koa,TSP,organic fraction'
              write(iout,*) Kp,eqkoa(i),TSP,fom
              call camxerr() 
            endif
         else                                           ! aerosol decrease
            daero = -min(abs(daero),ptconc(idxaer)-rtlbnd(nrtgas+i-1))
            if (abs(daero).ge.ptconc(idxaer)) then      ! Negative aerosol
              write(iout,*)'Error in PARTITIONRT:',
     &                     ' Partitioning more than available aerosol.'
              write(iout,*)'Available aerosol (ug/m3)',ptconc(idxaer)
              write(iout,*)'Mass subtracted (ug/m3)',abs(daero)
              write(iout,*)'Kp,Koa,TSP,organic fraction'
              write(iout,*) Kp,eqkoa(i),TSP,fom
              call camxerr() 
            endif
         endif
         ptconc(idx)    = ptconc(idx)    - daero/Mw     ! Update gas conc
         ptconc(idxaer) = ptconc(idxaer) + daero        ! Update aerosol conc
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
