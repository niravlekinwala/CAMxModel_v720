      subroutine hr_pan5(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_PAN5 solves PAN and RCO3 using Hertel's equations
c
c      Copyright 1996 - 2022
c     Ramboll
c     Created by the CMC version 5.2.6
c
c --- Subroutines Called:
c        none
c
c --- Called by:
c        EBISOLV
c
c --- Argument definitions:
c        y0  - initial Conc for this step (ppm)
c        yh  - current Conc iteration (ppm)
c        y1  - next Conc iteration (ppm)
c        H2O - water vapor Conc (ppm)
c        M   - total gas Conc (ppm)
c        O2  - oxygen Conc (ppm)
c        CH4 - methane Conc (ppm)
c        H2  - hydrogen Conc (ppm)
c        ny  - dimension of y0, y1 and yh
c        rk  - rate constants (ppm-n hr-1)
c        r   - reaction rates (hr-1)
c        nr  - dimension of rk and r
c        dt  - time step (hr)
c
c --- Includes:
      include "camx.prm"
      include "chmdat.inc"
      include "ddmchm.inc"
c
c --- Arguments:
      integer ny, nr
      real y0(ny+1), y1(ny+1), yh(ny+1)
      real rk(nr), r(nr)
      real H2O, M, O2, CH4, H2, dt
c
c --- Local variables:
      real N2
      real fwd, bck, self, A, B, C, Q
      real newRCO3, lsPAN, lsRCO3
c
c --- Entry Point
c
      N2 = M - O2
c
c
c --- New RCO3 excluding from PAN
c
      newRCO3 = dt*(0.0
     &       + ( 0.400)*r( 98)
     &       +          r( 99)
     &       +          r(101)
     &       +          r(102)
     &       +          r(103)
     &       +          r(104)
     &       +          r(106)
     &       +          r(107)
     &       + ( 2.000)*r(108)
     &       +          r(154)
     &       +          r(156)
     &       + ( 0.500)*r(157)
     &       + ( 0.500)*r(158)
     &       + ( 0.500)*r(159)
     &       +          r(160)
     &       +          r(161)
     &       +          r(162)
     &       +          r(163)
     &       +          r(208)
     &       +          r(210)
     &       + ( 0.620)*r(215)
     &       +          r(217)
     &       +          r(234)
     &       +          r(235)
     &       +          r(236)
     &       + ( 2.000)*r(237)
     &       + ( 0.305)*r(248)
     &       + ( 0.013)*r(253)
     &       + ( 0.340)*r(258)
     &       + ( 0.467)*r(266)
     &       + ( 0.400)*r(268)
     &       +          r(271)
     &       +          r(273)
     &       + ( 0.980)*r(279) )
c
c --- Loss of RCO3 excluding self-reaction
c
      lsRCO3 = 1.0 + dt*(0.0
     &       +          rk( 63)*yh(lNO2)
     &       +          rk( 66)*yh(lNO)
     &       +          rk( 67)*yh(lHO2)
     &       +          rk( 68)*yh(lNO3)
     &       +          rk( 69)*yh(lMEO2)
     &       +          rk( 70)*yh(lRO2C)
     &       +          rk( 71)*yh(lRO2X)
     &       +          rk( 82)*yh(lRCO3)
     &       +          rk( 93)*yh(lBZC3) )
c
c --- Loss of PAN
c
      lsPAN = 1.0 + dt*(0.0
     &       + rk( 64)
     &       + rk( 65) )
c
c
c --- Forward reaction of RCO3 to PAN
c
      fwd  = dt*(0.0
     &       +          rk( 63)*yh(lNO2) )
c
c --- Backward reactions of PAN to RCO3
c
      bck  = dt*(0.0
     &       +          rk( 64)
     &       + ( 0.600)*rk( 65) )
c
c --- RCO3 self-reaction
c
      self = dt*(0.0
     &       + ( 2.000)*rk( 72) )
c
c --- Solve for RCO3
c
      A = self*lsPAN
      B = ( lsPAN*lsRCO3) - (fwd*bck)
      C = lsPAN*( y0(lMCO3) + newRCO3 ) + bck*y0(lPAN)
      Q = -0.5 * ( B + SIGN(1.0,B)*SQRT(B*B + 4.0*A*C) )
c
c --- Update Concentrations
c
      y1(lMCO3) = MAX(1.0E-15, MAX(Q/A, -C/Q) )
      y1(lPAN)  = MAX(1.0E-15, (y0(lPAN) + fwd*y1(lMCO3)) / lsPAN )
c
      return
      end

