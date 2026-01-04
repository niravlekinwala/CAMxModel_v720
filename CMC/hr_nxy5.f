      subroutine hr_nxy5(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_NXY5 solves NO3 and N2O5 using Hertel's equations
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
      real self, fwd, bck
      real newNO3, lsNO3, lsN2O5, A, B, C, Q
c
c --- Entry Point
c
      N2 = M - O2
c
c
c --- New NO3 excluding from N2O5
c
      newNO3 = dt*(0.0
     &       +          r(  6)
     &       +          r(  8)
     &       +          r( 27)
     &       + ( 0.390)*r( 34)
     &       + ( 0.400)*r( 65)
     &       + ( 0.400)*r( 75)
     &       + ( 0.400)*r( 86)
     &       + ( 0.400)*r( 98) )
c
c --- Loss of NO3 excluding self-reaction(s)
c
      lsNO3 = 1.0 + dt*(0.0
     &       +          rk(  9)*yh(lNO)
     &       +          rk( 11)*yh(lNO2)
     &       +          rk( 15)*yh(lNO2)
     &       +          rk( 16)
     &       +          rk( 17)
     &       +          rk( 26)*yh(lOH)
     &       +          rk( 39)*yh(lHO2)
     &       +          rk( 49)*yh(lMEO2)
     &       +          rk( 54)*yh(lRO2C)
     &       +          rk( 59)*yh(lRO2X)
     &       +          rk( 68)*yh(lMCO3)
     &       +          rk( 78)*yh(lRCO3)
     &       +          rk( 89)*yh(lBZC3)
     &       +          rk(101)*yh(lMAC3)
     &       +          rk(207)*yh(lHCHO)
     &       +          rk(210)*yh(lCCHO)
     &       +          rk(213)*yh(lRCHO)
     &       +          rk(233)*yh(lGLY)
     &       +          rk(236)*yh(lMGLY)
     &       +          rk(239)*yh(lCRES)
     &       +          rk(245)*yh(lBALD)
     &       +          rk(256)*yh(lMACR)
     &       +          rk(265)*yh(lIPRD)
     &       +          rk(273)*yh(lGLYD)
     &       +          rk(276)*yh(lACRO)
     &       +          rk(516)*yh(lETHE)
     &       +          rk(520)*yh(lPRPE)
     &       +          rk(524)*yh(lBD13)
     &       +          rk(528)*yh(lISOP)
     &       +          rk(532)*yh(lAPIN)
     &       +          rk(550)*yh(lOLE1)
     &       +          rk(554)*yh(lOLE2)
     &       +          rk(560)*yh(lTERP)
     &       +          rk(564)*yh(lSESQ) )
c
c --- Loss of N2O5
c
      lsN2O5  = 1.0 + dt*(0.0
     &       +          rk( 12)
     &       +          rk( 13)*H2O
     &       +          rk( 14)*H2O*H2O )
c
c --- Forward reactions of NO3 to N2O5
c
      fwd  = dt*(0.0
     &       +          rk( 11)*yh(lNO2) )
c
c --- Backward reactions of N2O5 to NO3
c
      bck  = dt*(0.0
     &       +          rk( 12) )
c
c --- NO3 self-reaction excluding N2O5 production
c
      self = dt*(0.0
     &       + ( 2.000)*rk( 40) )
c
c --- Solve for NO3
c
      A = self*lsN2O5
      B = (lsN2O5*lsNO3) - (fwd*bck)
      C = lsN2O5*( y0(lNO3) + newNO3 ) + bck*y0(lN2O5)
      Q = -0.5 * ( B + SIGN(1.0,B)*SQRT(B*B + 4.0*A*C) )
c
c --- Update Concentrations
c
      y1(lNO3)  = MAX(1.0E-15, MAX(Q/A, -C/Q) )
      y1(lN2O5) = MAX(1.0E-15, (y0(lN2O5) + fwd*y1(lNO3)) / lsN2O5 )
c
      return
      end

