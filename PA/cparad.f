      subroutine cparad( rad, nr, pa, npa, nn, dtfact )
      use filunit
      use chmstry
      use camxcom
      use procan
      use tracer
c
c----CAMx v7.20 220430
c
c     CPARAD saves radical concentrations for CPA output
c     and sets the CPA names and position numbers on the 
c     first call from PASETUP
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        none
c
c     Input arguments:
c        rad                 radical concentrations array (ppm)
c        nr                  dimension of radical array
c        pa                  CPA parameter array (ppb units)
c        npa                 dimension of pa array
c        nn                  counter of CPA parameters filled
c        dtfact              ratio of time step to averaging interval
c
c     Output arguments:
c        pa                  CPA parameter array (ppb units)
c        nn                  counter of CPA parameters filled
c
c     Routines called:
c        none
c
c     Called by:
c        CHEMDRIV
c        PASETUP
c
      implicit none
      include 'camx.prm'
c
      integer npa, nr, nn, in_nn
      real    dtfact, ppbfact
      real    pa(npa), rad(nr)
c
c-----Entry point
c
      ppbfact = 1000.
      in_nn = nn
c
      nn = nn + 1
      ptname(nn)  = 'OH'
      cpadesc(nn)  = 'OH radical concentration'
      cpaunit(nn)  = 'hr^-1'
      pa(nn) = dtfact*rad(kOH)*ppbfact
c
      nn = nn + 1
      ptname(nn)  = 'HO2'
      cpadesc(nn)  = 'HO2 radical concentration'
      cpaunit(nn)  = 'hr^-1'
      pa(nn) = dtfact*rad(kHO2)*ppbfact
c
      return
      end
