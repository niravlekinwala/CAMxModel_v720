      subroutine pspgeo(iway,stdlon,cenlat,truelat1,xp,yp,lon,lat)
!
!----CAMx v7.20 220430
!
!     PSPGEO performs Polar Stereographic to geodetic (lat/lon) translation.
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
      real             :: stdlon   ! Longitude parallel to y-axis (-180->180E)
                                   ! Also used at (0,0) point
      real             :: cenlat   ! Latitude at (0,0) point (-90->+90)
      real             :: truelat1 ! True latitude (-90 -> 90 degrees N)
      real             :: xp       ! Cartesian X coordinate (km)
      real             :: yp       ! Cartesian Y coordinate (km)
      real             :: lat      ! Latitude (-90->90 deg N)
      real             :: lon      ! Longitude (-180->180 E)
!
!-----Locals
!
      real             :: hemi     ! 1 for NH, -1 for SH
      real             :: polex    ! Computed x-location of pole point
      real             :: poley    ! Computed y-location of pole point
      real             :: r0       ! Computed radius to (0,0) point
      real             :: rm
      real             :: ala,ala1
      real             :: alo,alo1
      real             :: reflon
      real             :: scale_top
      real             :: xx,yy
      real             :: gi2, r2
      real             :: arccos
!
!-----Entry point
!  
      if (truelat1 .LT. 0.) then
        hemi = -1.0
      else
        hemi = 1.0
      endif
!
!-----Compute the reference longitude by rotating 90 degrees to the east
!     to find the longitude line parallel to the positive x-axis.
!
      reflon = stdlon + 90.
! 
!-----Compute numerator term of map scale factor
!
      scale_top = 1. + hemi*SIN(truelat1*RAD_PER_DEG)
!  
!-----Compute radius to our known (0,0) point
!
      ala1 = cenlat*RAD_PER_DEG
      r0 = REARTH*COS(ala1)*scale_top/(1. + hemi*SIN(ala1))
      alo1 = (stdlon - reflon)*RAD_PER_DEG
      polex = -r0*COS(alo1)
      poley = -hemi*r0*SIN(alo1)
!
!-------------------------------------------------------------------------------
!     Calculate lat/lon from x/y
!-------------------------------------------------------------------------------
!
      if (iway .EQ. 1) then
!
!-----Compute radius to point of interest
!
        xx = xp - polex
        yy = (yp - poley)*hemi
        r2 = xx**2 + yy**2
! 
!-----Now the magic code
!
        if (r2 .EQ. 0.) then 
          lat = hemi*90.
          lon = reflon
        else
          gi2 = (REARTH*scale_top)**2.
          lat = DEG_PER_RAD*hemi*ASIN((gi2 - r2)/(gi2 + r2))
          arccos = ACOS(xx/SQRT(r2))
          if (yy .GT. 0) then
            lon = reflon + DEG_PER_RAD*arccos
          else
            lon = reflon - DEG_PER_RAD*arccos
          endif
        endif
!
!-----Convert to a -180 -> 180 East convention
!
        if (lon .GT.  180.) lon = lon - 360.
        if (lon .LT. -180.) lon = lon + 360.
!
!-------------------------------------------------------------------------------
!     Calculate x/y from lat/lon
!-------------------------------------------------------------------------------
!
      else
!
!-----Find radius to desired point
!
        ala = lat*RAD_PER_DEG
        rm = REARTH*COS(ala)*scale_top/(1. + hemi*SIN(ala))
        alo = (lon - reflon)*RAD_PER_DEG
        xp = polex + rm*COS(alo)
        yp = poley + hemi*rm*SIN(alo)
      endif
!
      return
      end 
