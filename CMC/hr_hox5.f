      subroutine hr_hox5(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_HOX5 solves the HOx family using Hertel's equations
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
      real t1, t2, t3, A, B, C, Q
      real newOH, newHO2, newHONO, newPNA
      real lsOH, lsHO2, lsHONO, lsPNA, lsO1D, yldOH, self
      real rOH_HO2, rHO2_OH, rOH_HONO, rHONO_OH, rHO2_PNA, rPNA_HO2
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
c --- Yield of OH from O1D
c
      yldOH  = ( 
     &         ( 2.000)*rk( 20)*H2O ) / lsO1D
c
c --- Conversion of HONO to OH
c
      rHONO_OH  = dt*(
     &       +          rk( 23) )
c
c --- Conversion of HO2 to OH
c
      rHO2_OH  = dt*(
     &       +          rk( 31)*yh(lNO)
     &       +          rk( 36)*yh(lO3)
     &       + ( 0.800)*rk( 39)*yh(lNO3) )
c
c --- Other OH production terms
c
      newOH = y0(lOH) + dt*(
     &       +          r( 28)
     &       + ( 0.390)*r( 34)
     &       + ( 2.000)*r( 41)
     &       +          r(124)
     &       +          r(126)
     &       + ( 0.500)*r(127)
     &       + ( 0.500)*r(128)
     &       + ( 0.500)*r(129)
     &       +          r(130)
     &       +          r(131)
     &       +          r(132)
     &       +          r(133)
     &       +          r(223)
     &       +          r(225)
     &       +          r(227)
     &       +          r(229)
     &       + ( 0.826)*r(247)
     &       + ( 0.826)*r(250)
     &       + ( 0.471)*r(253)
     &       + ( 0.208)*r(255)
     &       + ( 0.330)*r(258)
     &       + ( 0.164)*r(260)
     &       + ( 0.285)*r(264)
     &       + ( 0.330)*r(275)
     &       + ( 0.178)*r(278)
     &       +          r(280)
     &       +          r(282)
     &       + ( 0.160)*r(515)
     &       + ( 0.350)*r(519)
     &       + ( 0.080)*r(523)
     &       + ( 0.266)*r(527)
     &       + ( 0.728)*r(531)
     &       + ( 0.500)*r(535)
     &       + ( 0.128)*r(549)
     &       + ( 0.443)*r(553)
     &       + ( 0.499)*r(559)
     &       + ( 0.499)*r(563)
     &       +          r( 18)*yldOH )
c
c --- Conversion of PNA to HO2
c
      rPNA_HO2  = dt*(
     &       +          rk( 33)
     &       + ( 0.610)*rk( 34) )
c
c --- Conversion of OH to HO2
c
      rOH_HO2  = dt*(
     &       +          rk( 26)*yh(lNO3)
     &       +          rk( 29)*yh(lCO)
     &       +          rk( 30)*yh(lO3)
     &       +          rk( 42)*yh(lH2O2)
     &       +          rk( 44)*yh(lSO2)
     &       +          rk( 45)*H2
     &       +          rk(206)*yh(lHCHO)
     &       +          rk(218)*yh(lMEOH)
     &       +          rk(219)*yh(lFACD)
     &       + ( 0.148)*rk(228)*yh(lRAPX)
     &       + ( 0.630)*rk(232)*yh(lGLY)
     &       + ( 0.472)*rk(267)*yh(lPRD2)
     &       + ( 0.189)*rk(269)*yh(lRNO3)
     &       + ( 0.300)*rk(534)*yh(lACYE)
     &       + ( 0.570)*rk(536)*yh(lBENZ)
     &       + ( 0.181)*rk(537)*yh(lTOLU)
     &       + ( 0.159)*rk(538)*yh(lMXYL)
     &       + ( 0.161)*rk(539)*yh(lOXYL)
     &       + ( 0.159)*rk(540)*yh(lPXYL)
     &       + ( 0.022)*rk(541)*yh(lB124)
     &       + ( 0.950)*rk(542)*yh(lETOH)
     &       + ( 0.123)*rk(556)*yh(lARO1)
     &       + ( 0.077)*rk(557)*yh(lARO2) )
c
c --- Other HO2 production terms
c
      newHO2 = y0(lHO2) + dt*(
     &       +          r( 46)
     &       +          r( 49)
     &       + ( 2.000)*r( 51)
     &       + ( 0.500)*r( 55)
     &       + ( 0.500)*r( 60)
     &       + ( 0.900)*r( 69)
     &       +          r( 79)
     &       +          r( 90)
     &       +          r(102)
     &       +          r(114)
     &       +          r(116)
     &       + ( 0.500)*r(117)
     &       + ( 0.500)*r(118)
     &       + ( 0.500)*r(119)
     &       +          r(120)
     &       +          r(121)
     &       +          r(122)
     &       +          r(123)
     &       + ( 2.000)*r(204)
     &       +          r(207)
     &       +          r(209)
     &       +          r(212)
     &       +          r(223)
     &       +          r(225)
     &       + ( 0.142)*r(227)
     &       +          r(229)
     &       + ( 2.000)*r(230)
     &       + ( 0.630)*r(233)
     &       +          r(234)
     &       + ( 0.522)*r(247)
     &       + ( 1.023)*r(248)
     &       + ( 0.522)*r(250)
     &       + ( 0.554)*r(253)
     &       + ( 0.108)*r(255)
     &       + ( 0.670)*r(258)
     &       + ( 0.064)*r(260)
     &       + ( 0.400)*r(264)
     &       + ( 1.233)*r(266)
     &       + ( 0.344)*r(270)
     &       + ( 2.000)*r(272)
     &       + ( 0.830)*r(275)
     &       + ( 1.066)*r(278)
     &       +          r(485)
     &       + ( 0.500)*r(486)
     &       + ( 0.500)*r(487)
     &       + ( 0.500)*r(488)
     &       +          r(489)
     &       +          r(490)
     &       +          r(491)
     &       +          r(492)
     &       + ( 0.160)*r(515)
     &       + ( 0.800)*r(517)
     &       + ( 0.165)*r(519)
     &       + ( 0.080)*r(523)
     &       + ( 0.250)*r(525)
     &       + ( 0.066)*r(527)
     &       + ( 0.009)*r(531)
     &       + ( 1.500)*r(535)
     &       + ( 0.095)*r(549)
     &       + ( 0.094)*r(553)
     &       + ( 0.078)*r(559)
     &       + ( 0.078)*r(563) )
c
c
c --- Conversion of OH to HONO
c
      rOH_HONO  = dt*(
     &       +          rk( 22)*yh(lNO) )
c
c --- Other HONO production terms
c
      newHONO = y0(lHONO) + dt*(
     &       +          r(241) ) 
c
c --- Conversion of HO2 to PNA
c
      rHO2_PNA  = dt*(
     &       +          rk( 32)*yh(lNO2) )
c
c --- Other PNA production terms
c
      newPNA = y0(lPNA) 
c
c --- Net loss of OH
c
      lsOH = 1.0 + dt*(
     &       +          rk( 22)*yh(lNO)
     &       +          rk( 24)*yh(lHONO)
     &       +          rk( 25)*yh(lNO2)
     &       +          rk( 26)*yh(lNO3)
     &       +          rk( 27)*yh(lHNO3)
     &       +          rk( 29)*yh(lCO)
     &       +          rk( 30)*yh(lO3)
     &       +          rk( 35)*yh(lPNA)
     &       +          rk( 42)*yh(lH2O2)
     &       +          rk( 43)*yh(lHO2)
     &       +          rk( 44)*yh(lSO2)
     &       +          rk( 45)*H2
     &       +          rk(206)*yh(lHCHO)
     &       +          rk(208)*yh(lCCHO)
     &       +          rk(211)*yh(lRCHO)
     &       +          rk(214)*yh(lACET)
     &       +          rk(216)*yh(lMEK)
     &       +          rk(218)*yh(lMEOH)
     &       +          rk(219)*yh(lFACD)
     &       +          rk(220)*yh(lAACD)
     &       +          rk(221)*yh(lPACD)
     &       + ( 0.700)*rk(222)*yh(lCOOH)
     &       + ( 0.256)*rk(224)*yh(lROOH)
     &       + ( 0.160)*rk(226)*yh(lR6PX)
     &       + ( 0.861)*rk(228)*yh(lRAPX)
     &       +          rk(232)*yh(lGLY)
     &       +          rk(235)*yh(lMGLY)
     &       +          rk(238)*yh(lCRES)
     &       +          rk(240)*yh(lNPHE)
     &       +          rk(243)*yh(lBALD)
     &       +          rk(246)*yh(lAFG1)
     &       +          rk(249)*yh(lAFG2)
     &       +          rk(252)*yh(lAFG3)
     &       +          rk(254)*yh(lMACR)
     &       +          rk(259)*yh(lMVK)
     &       +          rk(263)*yh(lIPRD)
     &       +          rk(267)*yh(lPRD2)
     &       +          rk(269)*yh(lRNO3)
     &       +          rk(271)*yh(lGLYD)
     &       +          rk(274)*yh(lACRO)
     &       +          rk(279)*yh(lCO3H)
     &       +          rk(281)*yh(lRO3H)
     &       +          rk(513)*CH4
     &       +          rk(514)*yh(lETHE)
     &       +          rk(518)*yh(lPRPE)
     &       +          rk(522)*yh(lBD13)
     &       +          rk(526)*yh(lISOP)
     &       +          rk(530)*yh(lAPIN)
     &       + ( 0.300)*rk(534)*yh(lACYE)
     &       + ( 0.884)*rk(536)*yh(lBENZ)
     &       + ( 0.688)*rk(537)*yh(lTOLU)
     &       + ( 0.761)*rk(538)*yh(lMXYL)
     &       + ( 0.802)*rk(539)*yh(lOXYL)
     &       + ( 0.722)*rk(540)*yh(lPXYL)
     &       + ( 0.770)*rk(541)*yh(lB124)
     &       +          rk(542)*yh(lETOH)
     &       +          rk(543)*yh(lALK1)
     &       +          rk(544)*yh(lALK2)
     &       +          rk(545)*yh(lALK3)
     &       +          rk(546)*yh(lALK4)
     &       +          rk(547)*yh(lALK5)
     &       +          rk(548)*yh(lOLE1)
     &       +          rk(552)*yh(lOLE2)
     &       + ( 0.798)*rk(556)*yh(lARO1)
     &       + ( 0.822)*rk(557)*yh(lARO2)
     &       +          rk(558)*yh(lTERP)
     &       +          rk(562)*yh(lSESQ) )
c
c --- Loss of HO2, excluding self-reaction
c     (net either HO2 or OH produced)
c
      lsHO2 = 1.0 + rHO2_OH + rHO2_PNA + dt*(
     &       + ( 0.200)*rk( 39)*yh(lNO3)
     &       +          rk( 43)*yh(lOH)
     &       +          rk( 47)*yh(lMEO2)
     &       +          rk( 48)*yh(lMEO2)
     &       +          rk( 53)*yh(lRO2C)
     &       +          rk( 58)*yh(lRO2X)
     &       +          rk( 67)*yh(lMCO3)
     &       +          rk( 77)*yh(lRCO3)
     &       +          rk( 88)*yh(lBZC3)
     &       +          rk(100)*yh(lMAC3)
     &       +          rk(112)*yh(lBZO) )
c
c --- HO2 self-reaction
c
      self = dt*2.0*( 
     &       +          rk( 37)
     &       +          rk( 38)*H2O )
c
c --- Loss of HONO
c
      lsHONO = 1.0 + dt*(
     &       +          rk( 23)
     &       +          rk( 24)*yh(lOH) )
c
c --- Loss of PNA
c
      lsPNA = 1.0 + dt*(
     &       +          rk( 33)
     &       +          rk( 34)
     &       +          rk( 35)*yh(lOH) )
c
c --- Collect common terms
c
      t1 = 1.0 / ( lsOH*lsHONO - rHONO_OH*rOH_HONO )
      t2 = rOH_HO2*t1
      t3 = rPNA_HO2 / lsPNA
c
c --- Solve for HO2
c
      A = self
      B = lsHO2 - t3*rHO2_PNA - t2*rHO2_OH*lsHONO
      C = newHO2 + t3 * newPNA + t2*( newOH*lsHONO + newHONO*rHONO_OH )
      Q = -0.5 * (B + SIGN(1.0,B)*SQRT(B*B + 4.0*A*C))
c
c --- Update Concentrations
c
      y1(lHO2)  = MAX(1.0E-15, MAX(Q/A ,-C/Q) )
      y1(lOH)   = MAX(1.0E-15, ( ( newOH + rHO2_OH*y1(lHO2) )*lsHONO + 
     &                                        rHONO_OH*newHONO ) * t1 )
      y1(lPNA)  = MAX(1.0E-15, ( newPNA + rHO2_PNA*y1(lHO2) ) / lsPNA )
      y1(lHONO) = MAX(1.0E-15, ( newHONO + rOH_HONO*y1(lOH) ) / lsHONO )
c
      return
      end

