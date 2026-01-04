      subroutine getdelt(ncol,nrow,nlay,dx,dy,windu,windv,height,
     &                   mapscl,dtmx,kmax)
      use bndary
      use camxcom
c
c----CAMx v7.20 220430
c
c     GETDELT computes the time step size according to a maximum allowable
c     Courant number
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        11/30/99  Removed timestep calculation for horizontal diffusion
c        1/25/02   Made the CFL criterion a parameter
c        4/10/03   Now calculating max timestep by layer AND
c                  returning layer index of 'critical height'
c        1/16/06   Increased CFL to 0.9, reduced hmax to 10 m, revised kmax
c                  calculation to handle kmax = 1
c        5/25/06   CFL set according to advection solver
c        3/12/10   Increased Hmax for super-stepping to 2000 m (from 10 m)
c        11/3/15   Added map scale factor to timestep calculation
c        6/4/19    Using resultant X+Y winds instead of 1-D calculations
c        6/15/19   Revised application of map scale factor
c
c     Input arguments:
c        ncol                number of columne
c        nrow                number of rows
c        nlay                number of layers
c        dx                  cell sixe in x-direction (m)
c        dy                  cell size in y-direction (m)
c        windu               wind speed in x-direction (m/s)
c        windv               wind speed in y-direction (m/s)
c        height              layer height (m)
c        mapscl              map scale factor (unitless)
c
c     Output arguments:
c        dtmx                layer maximum time step size (s)
c        kmax                layer index of critical height hmax
c
c     Routines called:
c        none
c
c     Called by:
c        TIMESTEP
c
      implicit none
      include "camx.prm"
      include "flags.inc"
c
      integer ncol,nrow,nlay,kmax
      real dy,dx(nrow),windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &     height(ncol,nrow,nlay),mapscl(ncol,nrow)
      real dtmx(nlay)
c
      integer i,j,k,imid,jmid
      real hmax,cfl,tmpx,tmpy
      real uwind,vwind,uvwind
c
c-----Entry point
c
c
c  --- set height for super-stepping ---
c
      hmax = 9999999.
      if( lsuper ) hmax = 2000.
c
      cfl = 0.9
      if (iadvct.eq.3) cfl = 0.5
      do 30 k = 1,nlay
        dtmx(k) = 3600.
        tmpx = dtmx(k)
        tmpy = dtmx(k)
        do 20 j = 2,nrow-1
          do 10 i = 2,ncol-1
c
c-----Calculate maximum time step based on wind speed
c
            if (windu(i,j,k).gt.0. .and. windu(i-1,j,k).lt.0.) then
              uwind = abs(windu(i,j,k)) + abs(windu(i-1,j,k))
            else
              uwind = max(abs(windu(i,j,k)),abs(windu(i-1,j,k)))
            endif
            if (windv(i,j,k).gt.0. .and. windv(i,j-1,k).lt.0.) then
              vwind = abs(windv(i,j,k)) + abs(windv(i,j-1,k))
            else
              vwind = max(abs(windv(i,j,k)),abs(windv(i,j-1,k)))
            endif
            uvwind = sqrt(uwind**2 + vwind**2)
c           tmpx = cfl*dx(j)*mapscl(i,j)/(uvwind + 1.e-10)
c           tmpy = cfl*dy   *mapscl(i,j)/(uvwind + 1.e-10)
            tmpx = cfl*dx(j)/(uvwind + 1.e-10)/mapscl(i,j)
            tmpy = cfl*dy   /(uvwind + 1.e-10)/mapscl(i,j)
            dtmx(k) = min(dtmx(k),tmpx,tmpy)
 10       continue
 20     continue
 30   continue
c
c-----Determine layer index of 'critical height', which defines the depth
c     wherein the driving timestep for current grid will be determined 
c
      imid = ncol/2
      jmid = nrow/2
      do k = 2,nlay
        if (height(imid,jmid,k).gt.hmax) then 
          kmax = k-1
          kmax = max(1,kmax)
          goto 100
        endif
      enddo
      kmax = nlay
c
 100  return
      end
