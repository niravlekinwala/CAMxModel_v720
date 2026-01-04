      subroutine grdprep(nx,ny,cellon,cellat,mapscl,
     &                                  orgx,orgy,dx,dy,xsize,ysize)
      use grid
c 
c----CAMx v7.20 220430
c 
c     GRDPREP calculates coarse grid parameters for the following projections:
c        0. Lat/lon grid
c        1. Universal Transverse Mercator (UTM)
c        2. Lambert Conformal
c        3. Rotated (RAMS) Polar Stereographic
c        4. (WRF) Polar Stereographic 
c        5. (WRF) Mercator
c                           
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications: 
c        4/10/15      Added WRF Polar Stereographic and Mercator projections
c                     (Polar different from original Rotated Polar Stereographic
c                     from RAMS)
c 
c     Input arguments: 
c        nx                  number of columns
c        ny                  number of rows
c
c     Output arguments: 
c        cellon              array of cell centroid longitude (deg)
c        cellat              array of cell centroid latitude (deg)
c        mapscl              array of map-scale factors centered on each cell
c        orgx                origin of grid in X direction
c        orgy                origin of grid in X direction
c        dx                  cell width of grid in X direction
c        dy                  cell width of grid in Y direction
c        xsize               array of cell width in X direction
c        ysize               cell width in y direction
c             
c     Routines Called: 
c        UTMGEO
c        PSPGEO
c        LCPGEO
c             
c     Called by: 
c        STARTUP
c 
      implicit none
      include 'camx.prm'
      include 'flags.inc'
c
      integer nx,ny
      real cellon(nx,ny), cellat(nx,ny), mapscl(nx,ny)
      real orgx, orgy, dx, dy
      real xsize(ny), ysize
c
      integer i,j
      real deg2km,deg2rad
      real xloc,yloc,xloce,xlocw,alon,alat
      real colat,colat1,colat2,expn,xdel,ydel,elon,wlon,elat,wlat
c
      data deg2km /111.1338/
      data deg2rad /0.01745329/
c
c-----Entry point
c
c-----Cartesian coordinates are selected
c     Input cell sizes (dx and dy) are already in km
c
      if( .NOT. llatlon ) then
c
c-----Determine xsize and ysize (m) for coarse grid
c
        ysize = 1000.*dy
        do 10 j=1,ny
          xsize(j) = 1000.*dx
c
c-----Calculate lat/lon at cell centroids, and map-scale factor
c
          yloc = orgy + dy*(float(j) - 0.5)
          do 15 i = 1,nx
            xloc = orgx + dx*(float(i) - 0.5)
            if (lambrt) then
              call lcpgeo(1,polelon,polelat,tlat1,tlat2,xloc,yloc,alon,alat)
              colat1 = deg2rad*(90. - tlat1)
              colat = deg2rad*(90. - alat)
              if (tlat1 .ne. tlat2) then
                colat2 = deg2rad*(90. - tlat2)
                expn = (log(sin(colat1))    - log(sin(colat2))) /
     &                 (log(tan(colat1/2.)) - log(tan(colat2/2.)))
                mapscl(i,j) = sin(colat2)/sin(colat)*
     &                       (tan(colat/2.)/tan(colat2/2.))**expn
              else
                mapscl(i,j) = sin(colat1)/sin(colat)*
     &                        (tan(colat/2.)/tan(colat1/2.))**cos(colat1)
              endif
            elseif (lpolar) then
              call pspgeo(1,polelon,polelat,tlat1,xloc,yloc,alon,alat)
              mapscl(i,j) = (1. + sin(deg2rad*abs(tlat1)))/
     &                      (1. + sin(deg2rad*sign(1.,tlat1)*alat))
            elseif (lmerc) then
              call mrcgeo(1,polelon,polelat,tlat1,xloc,yloc,alon,alat)
              colat1 = deg2rad*(90. - tlat1)
              colat = deg2rad*(90. - alat)
              mapscl(i,j) = sin(colat1)/sin(colat)
            else
              if (lutm) then
                call utmgeo(1,iuzon,xloc,yloc,alon,alat)
                xloce = xloc + dx/2. 
                call utmgeo(1,iuzon,xloce,yloc,elon,elat)
                xlocw = xloc - dx/2. 
                call utmgeo(1,iuzon,xlocw,yloc,wlon,wlat)
              elseif (lrpolar) then
                call rpsgeo(1,polelon,polelat,xloc,yloc,alon,alat)
                xloce = xloc + dx/2. 
                call rpsgeo(1,polelon,polelat,xloce,yloc,elon,elat)
                xlocw = xloc - dx/2. 
                call rpsgeo(1,polelon,polelat,xlocw,yloc,wlon,wlat)
              endif
              xdel = deg2km*(elon - wlon)*cos(deg2rad*alat)  
              ydel = deg2km*(elat - wlat) 
              mapscl(i,j) = dx/sqrt(xdel*xdel + ydel*ydel) 
            endif
            cellon(i,j) = alon
            cellat(i,j) = alat
 15       continue
 10     continue
c
c-----Lat/lon coordinates are selected
c
      else
c
c-----Determine xsize and ysize (m) on coarse grid
c
        ysize = 1000.*deg2km*dy
        do 20 j = 1,ny
          alat = orgy + dy*(float(j) - 0.5)
          xsize(j) = 1000.*deg2km*dx*cos(alat*deg2rad)
c
c-----Calculate lat/lon at cell centroids
c
          do 25 i = 1,nx
            alon = orgx + dx*(float(i) - 0.5)
            cellon(i,j) = alon
            cellat(i,j) = alat
            mapscl(i,j) = 1.
 25        continue
 20     continue
      endif
c
      return
      end
