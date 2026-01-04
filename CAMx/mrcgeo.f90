      subroutine mrcgeo(iway,cenlon,cenlat,truelat1,xp,yp,lon,lat)
!
!----CAMx v7.20 220430
!
!     MRCGEO performs Mercator to geodetic (lat/lon) translation.
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
                                   ! 0 = geodetic to Polar Stereographic
                                   ! 1 = Polar Stereographic to geodetic
      real             :: cenlon   ! Longitude at (0,0) point (-180->180E)
      real             :: cenlat   ! Latitude at (0,0) point (-90->+90)
      real             :: truelat1 ! True latitude (-90 -> 90 degrees N)
      real             :: xp       ! Cartesian X coordinate (km)
      real             :: yp       ! Cartesian Y coordinate (km)
      real             :: lat      ! Latitude (-90->90 deg N)
      real             :: lon      ! Longitude (-180->180 E)
!
!-----Locals
!
      real             :: clain
      real             :: dlon
      real             :: deltalon
      real             :: r0
!
!-----Entry point
!
!-----Preliminary variables
!
      clain = COS(RAD_PER_DEG*truelat1)
      dlon = 1./(REARTH*clain)
!  
!-----Compute distance from equator to origin
!  
      r0 = 0.
      if (cenlat .NE. 0.) then
        r0 = (ALOG(TAN(0.5*((cenlat + 90.)*RAD_PER_DEG))))/dlon
      endif
!
!-------------------------------------------------------------------------------
!     Calculate lat/lon from x/y
!-------------------------------------------------------------------------------
!
      if (iway .EQ. 1) then
        lat = 2.*ATAN(EXP(dlon*(r0 + yp)))*DEG_PER_RAD - 90.
        lon = xp*dlon*DEG_PER_RAD + cenlon
        if (lon .GT.  180.) lon = lon - 360.
        if (lon .LT. -180.) lon = lon + 360.
!
!-------------------------------------------------------------------------------
!     Calculate x/y from lat/lon
!-------------------------------------------------------------------------------
!
      else
        deltalon = lon - cenlon
        if (deltalon .LT. -180.) deltalon = deltalon + 360.
        if (deltalon .GT.  180.) deltalon = deltalon - 360.
        xp = deltalon/(dlon*DEG_PER_RAD)
        yp = ALOG(TAN(0.5*((lat + 90.)*RAD_PER_DEG)))/dlon - r0
      endif
!
      return
      end
