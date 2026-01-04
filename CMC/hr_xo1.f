      subroutine hr_xo1(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_XO1 solves halogen oxides by EBI iterations
c
c      Copyright 1996 - 2022
c     Ramboll
c     Created by the CMC version 5.2.6 on 2020 02 20
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
      integer iter, nit
      real TINY
      parameter (TINY = 1.0E-15)
      real gainI2, lossI2
      real gainI, lossI
      real gainHOI, lossHOI
      real gainIO, lossIO
      real gainOIO, lossOIO
      real gainI2O2, lossI2O2
      real gainINO3, lossINO3
      real gainHIO3, lossHIO3
      real gainIXOY, lossIXOY
      real r209, r210, r211, r212, r213, r214, r215, r216, r217, r218
      real r219, r220, r221, r222, r223, r224
c
c --- Entry Point
c
      N2 = M - O2
c --- Perform iterations > 1
c     convergence will be tested in EBISOLV
c
      nit = 2
      if( y1(lIO) .GT. 250.0E-9) nit=3
c
      do iter = 2,nit
c
c --- Update reaction rates that change with iteration of y1
c
      r209 = rk(209)*y1(lI2)
      r210 = rk(210)*y1(lHOI)
      r211 = rk(211)*y1(lI)*yh(lO3)
      r212 = rk(212)*y1(lIO)
      r213 = rk(213)*y1(lIO)*y1(lIO)
      r214 = rk(214)*y1(lIO)*yh(lHO2)
      r215 = rk(215)*y1(lIO)*yh(lNO)
      r216 = rk(216)*y1(lIO)*yh(lNO2)
      r217 = rk(217)*y1(lOIO)
      r218 = rk(218)*y1(lOIO)*yh(lOH)
      r219 = rk(219)*y1(lOIO)*y1(lIO)
      r220 = rk(220)*y1(lOIO)*yh(lNO)
      r221 = rk(221)*y1(lI2O2)
      r222 = rk(222)*y1(lI2O2)*yh(lO3)
      r223 = rk(223)*y1(lINO3)
      r224 = rk(224)*y1(lINO3)*H2O
c
c --- Gain and loss terms for the XO species
c
      gainI2 = 0.0
c
      lossI2 = 0.0
     &       +          r209
c
      gainI = 0.0
     &       + ( 2.000)*r209
     &       +          r210
     &       +          r212
     &       + ( 0.400)*r213
     &       +          r215
     &       +          r217
     &       +          r221
     &       +          r223
c
      lossI = 0.0
     &       +          r211
c
      gainHOI = 0.0
     &       +          r214
     &       +          r224
c
      lossHOI = 0.0
     &       +          r210
c
      gainIO = 0.0
     &       +          r211
     &       +          r220
c
      lossIO = 0.0
     &       +          r212
     &       + ( 2.000)*r213
     &       +          r214
     &       +          r215
     &       +          r216
     &       +          r219
c
      gainOIO = 0.0
     &       + ( 0.400)*r213
     &       +          r221
c
      lossOIO = 0.0
     &       +          r217
     &       +          r218
     &       +          r219
     &       +          r220
c
      gainI2O2 = 0.0
     &       + ( 0.600)*r213
c
      lossI2O2 = 0.0
     &       +          r221
     &       +          r222
c
      gainINO3 = 0.0
     &       +          r216
c
      lossINO3 = 0.0
     &       +          r223
     &       +          r224
c
      gainHIO3 = 0.0
     &       +          r218
c
      lossHIO3 = 0.0
c
      gainIXOY = 0.0
     &       +          r219
     &       +          r222
c
      lossIXOY = 0.0
c
c --- EBI solution for the XO species
c
      y1(lI2) = MIN(1.0, MAX(TINY,(y0(lI2) + gainI2*dt)
     &                         / (1.0 + lossI2*dt/y1(lI2))))
c
      y1(lI) = MIN(1.0, MAX(TINY,(y0(lI) + gainI*dt)
     &                         / (1.0 + lossI*dt/y1(lI))))
c
      y1(lHOI) = MIN(1.0, MAX(TINY,(y0(lHOI) + gainHOI*dt)
     &                         / (1.0 + lossHOI*dt/y1(lHOI))))
c
      y1(lIO) = MIN(1.0, MAX(TINY,(y0(lIO) + gainIO*dt)
     &                         / (1.0 + lossIO*dt/y1(lIO))))
c
      y1(lOIO) = MIN(1.0, MAX(TINY,(y0(lOIO) + gainOIO*dt)
     &                         / (1.0 + lossOIO*dt/y1(lOIO))))
c
      y1(lI2O2) = MIN(1.0, MAX(TINY,(y0(lI2O2) + gainI2O2*dt)
     &                         / (1.0 + lossI2O2*dt/y1(lI2O2))))
c
      y1(lINO3) = MIN(1.0, MAX(TINY,(y0(lINO3) + gainINO3*dt)
     &                         / (1.0 + lossINO3*dt/y1(lINO3))))
c
      y1(lHIO3) = MIN(1.0, MAX(TINY,(y0(lHIO3) + gainHIO3*dt)
     &                         / (1.0 + lossHIO3*dt/y1(lHIO3))))
c
      y1(lIXOY) = MIN(1.0, MAX(TINY,(y0(lIXOY) + gainIXOY*dt)
     &                         / (1.0 + lossIXOY*dt/y1(lIXOY))))
c
      enddo
c
      return
      end

