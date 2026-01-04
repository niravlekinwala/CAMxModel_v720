c----CAMx v7.20 220430
c  
c     O3COLMAP.COM contains albedo and ozone column index maps and value classes
c                            
c      Copyright 1996 - 2022
c     Ramboll
c            
c     Modifications:  
c        02/11/11   Removed roughness panel
c        03/31/12   Removed haze and drought panels, haze classes to terrain
c                   height classes
c        08/08/14   Snow cover and albedo moved to camxfld.inc
c
c-----------------------------------------------------------------------
c     Variables to contain the albedo and ozone column data:
c
c     albcl   -- albedo column lookup table
c     ozcl    -- ozone column lookup table
c     trncl   -- terrain height lookup table
c     icdozn  -- index by cell/grid for ozone column
c     icdocn  -- optional index by cell/grid for ocean
c     lrdocn  -- flag by grid to show data read for ocean
c
c-----------------------------------------------------------------------
c
      real,    allocatable, dimension(:) :: albcl
      real,    allocatable, dimension(:) :: ozcl
      real,    allocatable, dimension(:) :: trncl
      integer, allocatable, dimension(:) :: icdozn
      integer, allocatable, dimension(:) :: icdocn
      logical, allocatable, dimension(:) :: lrdocn
c
