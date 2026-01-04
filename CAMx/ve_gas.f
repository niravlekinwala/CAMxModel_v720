      subroutine ve_gas(z0,deltaz,psih,ustar,diffrat,ve)
c
c----CAMx v7.20 220430
c
c     VE_GAS calculates an exchange/re-emission velocity for a specific gas 
c     species, grid cell, and land use category. 
c
c      Copyright 1996 - 2022
c     Ramboll
c 
c     Modifications: 
c        12/08/09   ---jjung---    Original development (for RTRAC)
c        04/29/13   --cemery---    Updated for species in the core model
c                                  Exchange across water-air interface ignored
c
c     Input arguments: 
c        z0                  surface roughness length (m)
c        deltaz              Layer 1 midpoint height (m)
c        psih                similarity stability correction term 
c        ustar               friction velocity (m/s)
c        diffrat             ratio of molecular diffusivity of water to species
c      
c     Output arguments: 
c        ve                  re-emission velocity (m/s)
c      
c     Routines called: 
c        none 
c      
c     Called by: 
c        SRFMOD
c 
      implicit none

      real z0,deltaz,psih,ustar,diffrat,ve
c
      real vk,rmin,d1,d2,vair,diffh2o
      real ra,rd,re
      real schmidt
c
      data vk/0.4/, rmin/1.0/, d1/2./, d2/0.667/
      data vair/1.5e-5/, diffh2o/2.3e-5/
c
c-----Entry point
c
      re = 0.
c
c-----Compute atmospheric resistance, RA
c
      ra = (alog(deltaz/z0) - psih)/(vk*ustar)
      ra = amax1(ra,rmin)
c
c-----Compute the deposition layer resistance, RD
c
      schmidt = vair*diffrat/diffh2o
      rd = d1*schmidt**d2/(vk*ustar)
      rd = amax1(rd,rmin)
c
c-----Re-emission velocity for this cell, land use, and species
c
 100  continue
      ve = 1./(ra + rd + re)
c
      return
      end
