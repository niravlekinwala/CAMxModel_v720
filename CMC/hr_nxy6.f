      subroutine hr_nxy6(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_NXY6 solves NO3 and N2O5 using Hertel's equations
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
     &       +          r(  5)
     &       +          r(  7)
     &       +          r( 29)
     &       + ( 0.390)*r( 51) )
c
c --- Loss of NO3 excluding self-reaction(s)
c
      lsNO3 = 1.0 + dt*(0.0
     &       +          rk( 14)
     &       +          rk( 15)
     &       +          rk( 16)*yh(lNO)
     &       +          rk( 17)*yh(lNO2)
     &       +          rk( 18)*yh(lNO2)
     &       +          rk( 46)*yh(lO)
     &       +          rk( 47)*yh(lOH)
     &       +          rk( 48)*yh(lHO2)
     &       +          rk( 49)*yh(lO3)
     &       +          rk( 78)*yh(lFORM)
     &       +          rk( 86)*yh(lALD2)
     &       +          rk(101)*yh(lALDX)
     &       +          rk(122)*yh(lOLE)
     &       +          rk(126)*yh(lETH)
     &       +          rk(130)*yh(lIOLE)
     &       +          rk(135)*yh(lCRES)
     &       +          rk(147)*yh(lISOP)
     &       +          rk(151)*yh(lISPD)
     &       +          rk(156)*yh(lTERP) )
c
c --- Loss of N2O5
c
      lsN2O5  = 1.0 + dt*(0.0
     &       +          rk( 19)*H2O
     &       +          rk( 20)*H2O*H2O
     &       +          rk( 21)
     &       +          rk( 53) )
c
c --- Forward reactions of NO3 to N2O5
c
      fwd  = dt*(0.0
     &       +          rk( 18)*yh(lNO2) )
c
c --- Backward reactions of N2O5 to NO3
c
      bck  = dt*(0.0
     &       +          rk( 21)
     &       +          rk( 53) )
c
c --- NO3 self-reaction excluding N2O5 production
c
      self = dt*(0.0
     &       + ( 2.000)*rk( 50) )
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

