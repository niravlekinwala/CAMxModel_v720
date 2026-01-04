      subroutine lcpgeo(iway,stdlon,cenlat,truelat1,truelat2,xp,yp,lon,lat)
!
!----CAMx v7.20 220430
!
!     LCPGEO performs Lambert Conformal to geodetic (lat/lon) translation.
!     Code based on MODULE_LLXY.F from WRF v3.3.
!
!      Portions Copyright 2022
!     Ramboll
!
!     Modifications:
!        none
!
!     Routines called:
!        none
!
!     Called by:
!        GRDPREP
!
      implicit none
      real, parameter :: PI = 3.141592653589793
      real, parameter :: DEG_PER_RAD = 180./PI
      real, parameter :: RAD_PER_DEG = PI/180.
      real, parameter :: REARTH = 6370.
!
!-----Arguments
!
      integer          :: iway     ! Conversion direction
                                   ! 0 = geodetic to Lambert Conformal
                                   ! 1 = Lambert Conformal to geodetic
      real             :: stdlon   ! Longitude parallel to y-axis (-180->180E)
                                   ! Also used at (0,0) point
      real             :: cenlat   ! Latitude at (0,0) point (-90->+90)
      real             :: truelat1 ! (-90 -> 90 degrees N)
      real             :: truelat2 !   "   "  "   "     "
      real             :: xp       ! Cartesian X coordinate (km)
      real             :: yp       ! Cartesian Y coordinate (km)
      real             :: lat      ! Latitude (-90->90 deg N)
      real             :: lon      ! Longitude (-180->180 E)
!
!-----Locals 
!
      real             :: polex    ! Computed x-location of pole point
      real             :: poley    ! Computed y-location of pole point
      real             :: cone     ! Cone factor for LC projections
      real             :: hemi     ! 1 for NH, -1 for SH
      real             :: r0       ! Computed radius to (0,0) point
      real             :: rr
      real             :: r2
      real             :: rm
      real             :: arg
      real             :: tl1r
      real             :: ctl1r
      real             :: xnew
      real             :: ynew
      real             :: chi,chi1,chi2
      real             :: xx
      real             :: yy
      real             :: deltalon
!
!-----Entry point
!
      if (truelat1 .LT. 0.) then
        hemi = -1.0
      else
        hemi = 1.0
      endif
      if (ABS(truelat2) .GT. 90.) then
        truelat2 = truelat1
      endif
!
!-----First, see if this is a secant or tangent projection.  For tangent
!     projections, truelat1 = truelat2 and the cone is tangent to the
!     Earth's surface at this latitude.  For secant projections, the cone
!     intersects the Earth's surface at each of the distinctly different
!     latitudes
!
      if (ABS(truelat1 - truelat2) .GT. 0.1) then
        cone = ALOG10(COS(truelat1*RAD_PER_DEG)) - &
               ALOG10(COS(truelat2*RAD_PER_DEG))
        cone = cone/(ALOG10(TAN((45. - ABS(truelat1)/2.)*RAD_PER_DEG)) - &
                     ALOG10(TAN((45. - ABS(truelat2)/2.)*RAD_PER_DEG)))
      else
         cone = SIN(ABS(truelat1)*RAD_PER_DEG)
      endif
!
!-----Convert truelat1 to radian and compute COS for later use
!
      tl1r = truelat1*RAD_PER_DEG
      ctl1r = COS(tl1r)
!
!-----Compute the radius to our known (0,0) point
!
      r0 = REARTH*ctl1r/cone * &
           (TAN((90.*hemi - cenlat)*RAD_PER_DEG/2.) / &
            TAN((90.*hemi - truelat1)*RAD_PER_DEG/2.))**cone
!     arg = cone*(deltalon1*RAD_PER_DEG)
!     polex = -hemi*r0*SIN(arg)
!     poley = r0*COS(arg)
      polex = 0.
      poley = r0
!
!-------------------------------------------------------------------------------
!     Calculate lat/lon from x/y
!-------------------------------------------------------------------------------
! 
      if (iway .EQ. 1) then
        chi1 = (90. - hemi*truelat1)*RAD_PER_DEG
        chi2 = (90. - hemi*truelat2)*RAD_PER_DEG
! 
!-----See if we are in the southern hemispere and flip the indices
!     if we are
!
        xnew = hemi*xp
        ynew = hemi*yp
! 
!-----Compute radius**2 to i/j location
!
        xx = xnew - polex
        yy = poley - ynew
        r2 = (xx*xx + yy*yy)
        rr = SQRT(r2)/REARTH
!
!-----Convert to lat/lon
!
        if (r2 .EQ. 0.) then
          lat = hemi*90.
          lon = stdlon
        else
          lon = stdlon + DEG_PER_RAD*ATAN2(hemi*xx,yy)/cone
          lon = MOD(lon + 360.,360.)
          if (lon .GT. +180.) lon = lon - 360.
          if (lon .LT. -180.) lon = lon + 360.
!
!-----Latitude determined by solving an equation adapted from:
!     Maling, D.H., 1973: Coordinate Systems and Map Projections
!     Equations #20 in Appendix I.  
!           
          if (chi1 .EQ. chi2) then
            chi = 2.*ATAN((rr/TAN(chi1))**(1./cone)*TAN(chi1*0.5))
          else
            chi = 2.*ATAN((rr*cone/SIN(chi1))**(1./cone)*TAN(chi1*0.5)) 
          endif
          lat = (90.0 - chi*DEG_PER_RAD)*hemi
        endif
!
!-------------------------------------------------------------------------------
!     Calculate x/y from lat/lon
!-------------------------------------------------------------------------------
! 
      else
!
!-----Compute deltalon between known longitude and standard lon and ensure
!     it is not in the cut zone
!
      deltalon = lon - stdlon
      if (deltalon .GT. +180.) deltalon = deltalon - 360.
      if (deltalon .LT. -180.) deltalon = deltalon + 360.
!      
!-----Radius to desired point
!
        rm = REARTH*ctl1r/cone * &
             (TAN((90.*hemi - lat)*RAD_PER_DEG/2.) / &
              TAN((90.*hemi - truelat1)*RAD_PER_DEG/2.))**cone
        arg = cone*(deltalon*RAD_PER_DEG)
        xp = polex + hemi*rm*SIN(arg)
        yp = poley - rm*COS(arg)
! 
!-----Finally, if we are in the southern hemisphere, flip the i/j
!     values to a coordinate system where (1,1) is the SW corner
!     (what we assume) which is different than the original NCEP
!     algorithms which used the NE corner as the origin in the 
!     southern hemisphere (left-hand vs. right-hand coordinate?)
!
        xp = hemi*xp
        yp = hemi*yp
      endif
!
      return
      end
