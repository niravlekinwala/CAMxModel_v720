      subroutine hr_nox5(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_NOX5 solves the NOx family using Hertel's equations
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
      real jO3_O1D, rNO2_O3P, rO3P_O3, rO3wNO
      real totO1D, totNO, totNO2, totO3P, totO3
      real lsO1D, lsNO, lsNO2, lsO3P, lsO3
      real yldO3P, rNO_NO2, rNO2_NO
      real t1, t2, t3, t4, t5, f1, f2, f3, f4, s1, s2, A, B, C, Q
c
c --- Entry Point
c
      N2 = M - O2
c
c --- O1D loss
c
      lsO1D = 
     &       + rk( 20)*H2O
     &       + rk( 21)*M
c
c --- Yield of O3P from O1D
c
      yldO3P  = ( rk( 21)*M ) / lsO1D
c
c --- Conversion of NO2 to NO
c     (rNO2_NO = rate of NO2 to NO)
c
      rNO2_NO  = dt*(
     &       +          rk(  1)
     &       +          rk(  5)*yh(lO3P) )
c
c --- NO production terms except rNO2_NO
c     (totNO = initial NO + new production)
c
      totNO = y0(lNO) + dt*(
     &       +          r( 15)
     &       +          r( 16)
     &       +          r( 23)
     &       + ( 0.000)*r(114)
     &       + ( 0.000)*r(124)
     &       + ( 0.000)*r(134)
     &       + ( 0.000)*r(144)
     &       + ( 0.000)*r(154)
     &       + ( 0.000)*r(164)
     &       + ( 0.000)*r(174)
     &       + ( 0.000)*r(184)
     &       + ( 0.000)*r(194)
     &       + ( 0.000)*r(283)
     &       + ( 0.000)*r(293)
     &       + ( 0.000)*r(303)
     &       + ( 0.000)*r(313)
     &       + ( 0.000)*r(323)
     &       + ( 0.000)*r(333)
     &       + ( 0.000)*r(343)
     &       + ( 0.000)*r(353)
     &       + ( 0.000)*r(363)
     &       + ( 0.000)*r(373)
     &       + ( 0.000)*r(383)
     &       + ( 0.000)*r(393)
     &       + ( 0.000)*r(403)
     &       + ( 0.000)*r(413)
     &       + ( 0.000)*r(423)
     &       + ( 0.000)*r(433)
     &       + ( 0.000)*r(443)
     &       + ( 0.000)*r(453)
     &       + ( 0.000)*r(463)
     &       + ( 0.000)*r(473)
     &       + ( 0.000)*r(483)
     &       + ( 0.000)*r(493)
     &       + ( 0.000)*r(503) )
c
c --- Net loss of NO
c     (net either NO or NO2 produced)
c
      lsNO = 1.0 + dt*(
     &       +          rk( 22)*yh(lOH)
     &       +          rk( 57)*yh(lRO2X) )
c
c --- Conversion of NO to NO2 except O3+NO (rO3wNO)
c
      rNO_NO2 = dt*(
     &       +          rk(  4)*yh(lO3P)
     &       +          rk(  9)*yh(lNO3)
     &       + ( 2.000)*rk( 10)*yh(lNO)*O2
     &       +          rk( 31)*yh(lHO2)
     &       +          rk( 46)*yh(lMEO2)
     &       +          rk( 52)*yh(lRO2C)
     &       +          rk( 66)*yh(lMCO3)
     &       +          rk( 76)*yh(lRCO3)
     &       +          rk( 87)*yh(lBZC3)
     &       +          rk( 99)*yh(lMAC3)
     &       +          rk(134)*yh(lXNO2) )
c
c --- Remaining NO2 production
c
      totNO2 = y0(lNO2) + dt*(
     &       +          r(  9)
     &       +          r( 12)
     &       +          r( 17)
     &       +          r( 24)
     &       +          r( 26)
     &       +          r( 28)
     &       +          r( 33)
     &       + ( 0.610)*r( 34)
     &       +          r( 35)
     &       + ( 0.800)*r( 39)
     &       + ( 2.000)*r( 40)
     &       +          r( 49)
     &       +          r( 54)
     &       +          r( 59)
     &       +          r( 64)
     &       + ( 0.600)*r( 65)
     &       +          r( 68)
     &       +          r( 74)
     &       + ( 0.600)*r( 75)
     &       +          r( 78)
     &       +          r( 85)
     &       + ( 0.600)*r( 86)
     &       +          r( 89)
     &       +          r( 97)
     &       + ( 0.600)*r( 98)
     &       +          r(101)
     &       +          r(136)
     &       + ( 0.500)*r(137)
     &       + ( 0.500)*r(138)
     &       + ( 0.500)*r(139)
     &       +          r(140)
     &       +          r(141)
     &       +          r(142)
     &       +          r(143)
     &       + ( 0.019)*r(269)
     &       +          r(270) )
c
c --- Net loss of NO2
c     (net either NO or NO2 produced)
c
      lsNO2 = 1.0 + dt*(
     &       +          rk(  6)*yh(lO3P)
     &       +          rk(  8)*yh(lO3)
     &       +          rk( 11)*yh(lNO3)
     &       +          rk( 25)*yh(lOH)
     &       +          rk( 32)*yh(lHO2)
     &       +          rk( 63)*yh(lMCO3)
     &       +          rk( 73)*yh(lRCO3)
     &       +          rk( 84)*yh(lBZC3)
     &       +          rk( 96)*yh(lMAC3)
     &       +          rk(109)*yh(lTBUO)
     &       +          rk(111)*yh(lBZO) )
c
c --- Production of O3P except NO2+hv
c
      totO3P = y0(lO) + dt*(
     &       +          r( 17)
     &       +          r( 18)*yldO3P
     &       +          r( 19) )
c
c --- Net loss of O3P
c
      lsO3P = 1.0 + dt*(
     &       +          rk(  2)*O2*M
     &       +          rk(  3)*yh(lO3)
     &       +          rk(  4)*yh(lNO)
     &       +          rk(  5)*yh(lNO2)
     &       +          rk(  6)*yh(lNO2)
     &       +          rk(257)*yh(lMACR)
     &       +          rk(261)*yh(lMVK)
     &       +          rk(277)*yh(lACRO)
     &       +          rk(517)*yh(lETHE)
     &       +          rk(521)*yh(lPRPE)
     &       +          rk(525)*yh(lBD13)
     &       +          rk(529)*yh(lISOP)
     &       +          rk(533)*yh(lAPIN)
     &       +          rk(551)*yh(lOLE1)
     &       +          rk(555)*yh(lOLE2)
     &       +          rk(561)*yh(lTERP)
     &       +          rk(565)*yh(lSESQ) )
c
c --- Production of O3 except O3P+O2
c
      totO3 = y0(lO3) + dt*(
     &       + ( 0.300)*r( 67)
     &       + ( 0.250)*r( 77)
     &       + ( 0.250)*r( 88)
     &       + ( 0.250)*r(100) )
c
c --- Net loss of O3 except O3+NO (rO3wNO)
c
      lsO3 = 1.0 + dt*(
     &       +          rk(  3)*yh(lO3P)
     &       +          rk(  8)*yh(lNO2)
     &       +          rk( 18)
     &       +          rk( 19)
     &       +          rk( 30)*yh(lOH)
     &       +          rk( 36)*yh(lHO2)
     &       +          rk(247)*yh(lAFG1)
     &       +          rk(250)*yh(lAFG2)
     &       +          rk(253)*yh(lAFG3)
     &       +          rk(255)*yh(lMACR)
     &       +          rk(260)*yh(lMVK)
     &       +          rk(264)*yh(lIPRD)
     &       +          rk(275)*yh(lACRO)
     &       +          rk(515)*yh(lETHE)
     &       +          rk(519)*yh(lPRPE)
     &       +          rk(523)*yh(lBD13)
     &       +          rk(527)*yh(lISOP)
     &       +          rk(531)*yh(lAPIN)
     &       +          rk(535)*yh(lACYE)
     &       +          rk(549)*yh(lOLE1)
     &       +          rk(553)*yh(lOLE2)
     &       +          rk(559)*yh(lTERP)
     &       +          rk(563)*yh(lSESQ) )
c
c --- Specific reactions
c
      jO3_O1D  = rk( 18)
      rNO2_O3P = dt*rk(  1)
      rO3P_O3  = dt*rk(  2)*O2*M
      rO3wNO   = dt*rk(  7)
c
c --- Collect common terms
c
      t1 = rNO2_O3P / lsNO2
      t2 = rNO2_NO / lsNO2
      t3 = rNO_NO2 / lsNO
      t4 = rO3P_O3  / lsO3P
      t5 = t3*totNO - t2*totNO2
      f1 = 1.0 + t2 + t3
      f2 = t1*t4
      f3 = lsO3*lsNO + rO3wNO*totNO
      f4 = totO3 + totO3P * t4
c
c --- Solve for change in NO and NO2 
c
      A = rO3wNO * (f1 - f2)
      B = f1*f3 + rO3wNO*(f2*(totNO2 - totNO) + f4 + t5)
      C = rO3wNO*totNO*(f4 + totNO2 * f2) + f3*t5
      Q = -0.5 * (B + SIGN(1.0,B)*SQRT(B*B - 4.0*A*C))
      Q = MAX(Q/A ,C/Q)
c
c --- Update Concentrations
c
      y1(lNO)  = MAX(1.0E-15, (totNO + Q) / lsNO )
      y1(lNO2) = MAX(1.0E-15, (totNO2 - Q) / lsNO2 )
      s1 = totO3P + rNO2_O3P*y1(lNO2)
      s2 = t4*s1
      y1(lO3)  = MAX(1.0E-15, (totO3 + s2) / 
     &                             ( lsO3 + rO3wNO * y1(lNO) ) )
      y1(lO3P)   = MAX(1.0E-15, s1 / lsO3P )
      totO1D = jO3_O1D*y1(lO3)
      y1(lO1D) = MAX(1.0E-15, totO1D / lsO1D )
c
      return
      end

