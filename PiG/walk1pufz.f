      subroutine walk1pufz(dt,nx,ny,nz,ii,jj,zpig,height,ww)
c     use grid
c
c----CAMx v7.20 220430
c
c     WALK1PUFZ calculates a new puff vertical location by integrating
c     dz/w(z) = dt through a vertical column, where
c     w(z) = w(k-1) + (w(k) - w(k-1))*(z/depth).
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        10/19/16        Bug fix for puffs reaching layer interface when w=0
c
c     Input arguments:
c        dt                  remaining timestep (s)
c        nx                  number of cells in x-direction
c        ny                  number of cells in y-direction
c        nz                  number of layers
c        ii                  i-cell index of puff location 
c        jj                  j-cell index of puff location 
c        zpig                z-coord of puff location (m)
c        height              gridded layer height (m)
c        ww                  vertical velocity profile (m/s)
c
c     Output arguments:
c        zpig                z-coord of puff location (m)
c
c     Subroutines Called:
c        NONE
c
c     Called by:
c        PIGWALK
c
      implicit none
      include "camx.prm"
c
      integer nx,ny,nz,ii,jj
      real zpig,dt
      real windu(nx,ny,nz),windv(nx,ny,nz),height(nx,ny,nz),ww(0:nz)
      real depth,zcel0,zcell,az,bz,dtz,dttmp
      integer kount,iedge,kp,kk
      real zfrac
c
c-----Entry point
c
      dttmp = dt
      kount = 0
 10   iedge = 0
      kount = kount + 1
c
c-----Calculate distance between puff midpoint and layer interface below
c
      kp = nz
      do kk = 1,nz
        if (height(ii,jj,kk) .GE. zpig) then
          kp = kk
          goto 11
        endif
      enddo

 11   continue
      depth = height(ii,jj,kp)
      zcel0 = zpig
      if (kp.gt.1) then
        depth = height(ii,jj,kp) - height(ii,jj,kp-1)
        zcel0 = zpig - height(ii,jj,kp-1)
      endif
c
c-----Compute the new position after the whole time step
c
      iedge = 0
      az = ww(kp-1)
      bz = (ww(kp) - ww(kp-1))/depth
      if (abs(bz).lt.1.e-6) then
        zcell = zcel0 + az*dttmp
      else
        zcell = (-az + (az + bz*zcel0)*exp(bz*dttmp))/bz
      endif
c
c-----If the new position is beyond the layer boundaries, compute
c     the time needed to reach the boundary
c
      dtz = dttmp
      if (zcell.gt.depth) then
        iedge = 1
        if (abs(bz).lt.1.e-6) then
          dtz = (depth*1.001 - zcel0)/az
        else
          zfrac = (az + bz*depth*1.001)/(az + bz*zcel0)
          if (zfrac.le.0.) then
            iedge = 0
          else
            dtz = alog(zfrac)/bz
            if (dtz.gt.dttmp) then
              iedge = 0
              dtz = dttmp
            endif
          endif
        endif
      elseif (zcell.lt.0.) then
        iedge = 1
        if (abs(bz).lt.1.e-6) then
          dtz = (-depth*0.001 - zcel0)/az
        else
          zfrac = (az - bz*depth*0.001)/(az + bz*zcel0)
          if (zfrac.le.0.) then
            iedge = 0
          else
            dtz = alog(zfrac)/bz
            if (dtz.gt.dttmp) then
              iedge = 0
              dtz = dttmp
            endif
          endif
        endif
      endif
c
c-----If the puff reaches a boundary, compute the remaining time and new
c     position
c
      if (iedge.eq.1) then
        if (abs(bz).lt.1.e-6) then
          zcell = zcel0 + az*dtz
        else
          zcell = (-az + (az + bz*zcel0)*exp(bz*dtz))/bz
        endif
        dttmp = dttmp - dtz
      else
        dttmp = 0.
      endif
c
      zpig = max(0.,zcell)
      if (kp.gt.1) zpig = min(height(ii,jj,nz),
     &                        zcell + height(ii,jj,kp-1))

      if (kount.le.50 .and. dttmp.gt.0.) goto 10
c
      return
      end
