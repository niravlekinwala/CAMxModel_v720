      subroutine hr_pan4(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_PAN4 solves PAN and RCO3 using Hertel's equations
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
     &       +          r( 95)
     &       +          r(104)
     &       +          r(105)
     &       + ( 0.800)*r(110)
     &       +          r(112)
     &       +          r(116)
     &       +          r(117)
     &       +          r(118)
     &       + ( 0.500)*r(125)
     &       + ( 0.620)*r(126)
     &       +          r(127)
     &       + ( 0.269)*r(151)
     &       + ( 0.143)*r(152)
     &       + ( 0.208)*r(154)
     &       + ( 0.300)*r(189)
     &       + ( 0.600)*r(191)
     &       + ( 0.120)*r(195) )
c
c --- Loss of RCO3 excluding self-reaction
c
      lsRCO3 = 1.0 + dt*(0.0
     &       +          rk( 53)*yh(lNO)
     &       +          rk( 54)*yh(lNO2)
     &       +          rk( 57)*yh(lHO2)
     &       +          rk( 60)*yh(lCXO3)
     &       +          rk( 73)*yh(lMEO2)
     &       +          rk( 77)*yh(lXO2H)
     &       +          rk( 81)*yh(lXO2)
     &       +          rk( 85)*yh(lXO2N)
     &       +          rk(146)*yh(lISO2)
     &       +          rk(161)*yh(lEPX2)
     &       +          rk(169)*yh(lBZO2)
     &       +          rk(174)*yh(lTO2)
     &       +          rk(180)*yh(lXLO2)
     &       +          rk(203)*yh(lOPO3) )
c
c --- Loss of PAN
c
      lsPAN = 1.0 + dt*(0.0
     &       + rk( 55)
     &       + rk( 56) )
c
c
c --- Forward reaction of RCO3 to PAN
c
      fwd  = dt*(0.0
     &       +          rk( 54)*yh(lNO2) )
c
c --- Backward reactions of PAN to RCO3
c
      bck  = dt*(0.0
     &       +          rk( 55)
     &       + ( 0.600)*rk( 56) )
c
c --- RCO3 self-reaction
c
      self = dt*(0.0
     &       + ( 2.000)*rk( 59) )
c
c --- Solve for RCO3
c
      A = self*lsPAN
      B = ( lsPAN*lsRCO3) - (fwd*bck)
      C = lsPAN*( y0(lC2O3) + newRCO3 ) + bck*y0(lPAN)
      Q = -0.5 * ( B + SIGN(1.0,B)*SQRT(B*B + 4.0*A*C) )
c
c --- Update Concentrations
c
      y1(lC2O3) = MAX(1.0E-15, MAX(Q/A, -C/Q) )
      y1(lPAN)  = MAX(1.0E-15, (y0(lPAN) + fwd*y1(lC2O3)) / lsPAN )
c
      return
      end

