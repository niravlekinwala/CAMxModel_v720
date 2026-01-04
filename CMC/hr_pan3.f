      subroutine hr_pan3(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_PAN3 solves PAN and RCO3 using Hertel's equations
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
     &       +          r(105)
     &       +          r(106)
     &       +          r(107)
     &       + ( 0.800)*r(113)
     &       +          r(115)
     &       +          r(119)
     &       +          r(120)
     &       +          r(121)
     &       + ( 0.500)*r(128)
     &       + ( 0.620)*r(129)
     &       +          r(130)
     &       + ( 0.269)*r(158)
     &       + ( 0.143)*r(159)
     &       + ( 0.208)*r(161)
     &       + ( 0.300)*r(197)
     &       + ( 0.600)*r(199)
     &       + ( 0.120)*r(203)
     &       +          r(261) )
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
     &       +          rk(153)*yh(lISO2)
     &       +          rk(168)*yh(lEPX2)
     &       +          rk(177)*yh(lBZO2)
     &       +          rk(182)*yh(lTO2)
     &       +          rk(188)*yh(lXLO2)
     &       +          rk(211)*yh(lOPO3) )
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

