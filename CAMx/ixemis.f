      subroutine ixemis(dx,dy,mapscl,iocn,fwtr,o3,tsrf,wind,i2flx,
     &                  hoiflx)
c
c     use chmstry
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    This routine calculates the amount of inorganic iodine source gas
c    emissions (I2 and HOI) from ocean surfaces as a function of O3,
c    surface temperature and wind speed.
c    Equations from Prados-Roman et al. (2015), Atm. Chem. Phys.,
c                   doi:10.5194/acp-15-2215-2015.
c
c    Copyright 2022
c    Ramboll
c
c    Input arguments:
c     dx        R  E-W grid cell size (m)
c     dy        R  N-S grid cell size (m)
c     mapscl    R  map scale factor
c     iocn      I  Ocean flag for this grid cell
c     fwtr      R  Cell water fraction
c     o3        R  Ozone concentration (ppb)
c     tsrf      R  Surface temperature (K)
c     wind      R  Layer 1 wind speed (m/s)
c
c   Output arguments:
c     i2flx     R  I2 oceanic emissions flux (mol/s)
c     hoiflx    R  HOI oceanic emissions flux (mol/s)
c
c   Routines called:
c       None
c
c   Called by:
c       EMISS 
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     12/15/15   ---cemery---    Original development
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      implicit none
c     include 'camx.prm'
c     include 'deposit.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer iocn
      real    dx, dy
      real    mapscl
      real    o3
      real    fwtr
      real    tsrf
      real    wind
      real    i2flx
      real    hoiflx
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      real      iaq 
      real      area
      real      nano
      real      secday
c
      data nano /1.e-9/
      data secday /86400./
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c-----Calculate emission flux of I2 and HOI for this cell
c
      i2flx = 0.
      hoiflx = 0.
      if (fwtr.GT.0. .AND. iocn.GT.0) then
        iaq = 1.46e6*EXP(-9134./tsrf)
        i2flx = o3*(iaq**1.3)*(1.74e9 - 6.54e8*alog(wind)) ! (nmol/m2/d)
        hoiflx = o3*(4.15e5*sqrt(iaq)/wind - 20.6/wind - 
     &               2.36e4*sqrt(iaq))                     ! (nmol/m2/d)
c
c-----Units conversion, cell/ocean fraction, lower bound is zero
c
        area = fwtr*dx*dy/mapscl**2.
        i2flx = max(0.,i2flx)*nano*area/secday
        hoiflx = max(0.,hoiflx)*nano*area/secday           ! (mol/s)
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
