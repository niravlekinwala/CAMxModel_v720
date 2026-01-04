      subroutine hr_hox7(y0,y1,yh,H2O,M,O2,CH4,H2,ny,rk,r,nr,dt)
      implicit none
c
c----CAMx v7.20 220430
c
c     HR_HOX7 solves the HOx family using Hertel's equations
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
     &       + rk( 10)*M
     &       + rk( 11)*H2O
c
c --- Yield of OH from O1D
c
      yldOH  = ( 
     &         ( 2.000)*rk( 11)*H2O ) / lsO1D
c
c --- Conversion of HONO to OH
c
      rHONO_OH  = dt*(
     &       +          rk( 39) )
c
c --- Conversion of HO2 to OH
c
      rHO2_OH  = dt*(
     &       +          rk( 13)*yh(lO3)
     &       +          rk( 15)*yh(lO)
     &       +          rk( 25)*yh(lNO)
     &       +          rk( 32)*yh(lNO3)
     &       + ( 0.500)*rk( 60)*yh(lC2O3)
     &       + ( 0.250)*rk( 67)*yh(lCXO3)
     &       + ( 0.500)*rk( 73)*yh(lOPO3)
     &       + ( 0.060)*rk(180)*yh(lISO2)
     &       + ( 1.700)*rk(193)*yh(lEPX2)
     &       + ( 0.350)*rk(198)*yh(lAPO2)
     &       + ( 0.070)*rk(204)*yh(lTPO2) )
cgwilson-CB7     &       + ( 0.100)*rk(210)*yh(lSQO2) )
c
c --- Other OH production terms
c
      newOH = y0(lOH) + dt*(
     &       + ( 2.000)*r( 21)
     &       +          r( 23)
     &       +          r( 44)
     &       + ( 0.410)*r( 47)
     &       +          r( 83)
     &       +          r( 94)
     &       + ( 0.190)*r(111)
     &       + ( 0.170)*r(140)
     &       + ( 0.334)*r(143)
     &       + ( 0.500)*r(146)
     &       + ( 0.500)*r(162)
     &       + ( 0.500)*r(166)
     &       + ( 1.700)*r(182)
     &       + ( 0.280)*r(183)
     &       + ( 0.250)*r(184)
     &       + ( 0.200)*r(187)
     &       + ( 1.650)*r(189)
     &       + ( 0.700)*r(192)
     &       + ( 0.850)*r(194)
     &       + ( 0.770)*r(200)
     &       + ( 0.700)*r(201)
     &       + ( 0.470)*r(206)
     &       + ( 0.480)*r(207)
     &       + ( 0.080)*r(212)
     &       + ( 0.470)*r(213)
     &       +          r(215)
     &       +          r(  9)*yldOH )
c
c --- Conversion of PNA to HO2
c
      rPNA_HO2  = dt*(
     &       +          rk( 46)
     &       + ( 0.590)*rk( 47) )
c
c --- Conversion of OH to HO2
c
      rOH_HO2  = dt*(
     &       +          rk( 12)*yh(lO3)
     &       +          rk( 14)*yh(lO)
     &       +          rk( 22)*yh(lH2O2)
     &       +          rk( 31)*yh(lNO3)
     &       +          rk( 49)*H2
     &       +          rk( 50)*yh(lCO)
     &       +          rk( 51)*yh(lSO2)
     &       +          rk( 98)*yh(lMEOH)
     &       + ( 0.900)*rk( 99)*yh(lETOH)
     &       +          rk(100)*yh(lFORM)
     &       + ( 0.200)*rk(110)*yh(lGLYD)
     &       +          rk(113)*yh(lGLY)
     &       + ( 0.610)*rk(122)*yh(lKET)
     &       +          rk(123)*yh(lHACT)
     &       +          rk(124)*yh(lFACD)
     &       + ( 0.300)*rk(138)*yh(lETHY)
     &       + ( 0.530)*rk(148)*yh(lBENZ)
     &       + ( 0.180)*rk(152)*yh(lTOL)
     &       + ( 0.155)*rk(156)*yh(lXYL)
     &       +          rk(168)*yh(lCRES)
     &       + ( 0.200)*rk(175)*yh(lCAT1)
     &       + ( 0.100)*rk(185)*yh(lISPD)
     &       + ( 0.500)*rk(188)*yh(lISPX)
     &       + ( 0.200)*rk(191)*yh(lEPOX) )
c
c --- Other HO2 production terms
c
      newHO2 = y0(lHO2) + dt*(
     &       +          r( 23)
     &       + ( 0.350)*r( 68)
     &       + ( 0.800)*r( 69)
     &       +          r( 78)
     &       + ( 0.900)*r( 80)
     &       + ( 0.370)*r( 81)
     &       +          r( 84)
     &       + ( 0.600)*r( 86)
     &       +          r( 94)
     &       + ( 2.000)*r(101)
     &       +          r(103)
     &       +          r(106)
     &       +          r(109)
     &       + ( 1.400)*r(111)
     &       + ( 2.000)*r(114)
     &       +          r(115)
     &       +          r(116)
     &       + ( 0.960)*r(121)
     &       + ( 0.620)*r(136)
     &       +          r(137)
     &       + ( 0.270)*r(140)
     &       + ( 0.080)*r(143)
     &       + ( 0.918)*r(149)
     &       +          r(151)
     &       + ( 0.860)*r(153)
     &       +          r(155)
     &       + ( 0.860)*r(157)
     &       +          r(159)
     &       +          r(160)
     &       + ( 0.560)*r(162)
     &       + ( 0.700)*r(164)
     &       +          r(174)
     &       + ( 0.900)*r(179)
     &       + ( 0.350)*r(182)
     &       + ( 0.500)*r(183)
     &       + ( 0.200)*r(189)
     &       +          r(192)
     &       +          r(194)
     &       + ( 0.770)*r(197)
     &       + ( 0.500)*r(199)
     &       + ( 0.170)*r(200)
     &       + ( 0.700)*r(201)
     &       + ( 0.750)*r(203)
     &       + ( 0.500)*r(205)
     &       + ( 0.160)*r(206)
     &       + ( 0.300)*r(207)
     &       + ( 0.600)*r(209)
     &       + ( 0.500)*r(211)
     &       + ( 0.080)*r(212) )
c
c
c --- Conversion of OH to HONO
c
      rOH_HONO  = dt*(
     &       +          rk( 38)*yh(lNO) )
c
c --- Other HONO production terms
c
      newHONO = y0(lHONO) + dt*(
     &       +          r(174) ) 
c
c --- Conversion of HO2 to PNA
c
      rHO2_PNA  = dt*(
     &       +          rk( 45)*yh(lNO2) )
c
c --- Other PNA production terms
c
      newPNA = y0(lPNA) 
c
c --- Net loss of OH
c
      lsOH = 1.0 + dt*(
     &       +          rk( 12)*yh(lO3)
     &       +          rk( 14)*yh(lO)
     &       +          rk( 16)*yh(lOH)
     &       +          rk( 17)*yh(lOH)
     &       +          rk( 18)*yh(lHO2)
     &       +          rk( 22)*yh(lH2O2)
     &       +          rk( 31)*yh(lNO3)
     &       +          rk( 38)*yh(lNO)
     &       +          rk( 40)*yh(lHONO)
     &       +          rk( 41)*yh(lNO2)
     &       +          rk( 42)*yh(lNO2)*H2O
     &       +          rk( 43)*yh(lHNO3)
     &       +          rk( 48)*yh(lPNA)
     &       +          rk( 49)*H2
     &       +          rk( 50)*yh(lCO)
     &       +          rk( 51)*yh(lSO2)
     &       +          rk( 53)*yh(lDMS)
     &       +          rk( 54)*yh(lDMS)*O2
     &       +          rk( 66)*yh(lPANX)
     &       +          rk( 72)*yh(lOPAN)
     &       + ( 0.600)*rk( 82)*yh(lMEPX)
     &       + ( 0.600)*rk( 93)*yh(lROOH)
     &       +          rk( 95)*yh(lNTR1)
     &       +          rk( 98)*yh(lMEOH)
     &       +          rk( 99)*yh(lETOH)
     &       +          rk(100)*yh(lFORM)
     &       +          rk(104)*yh(lALD2)
     &       +          rk(107)*yh(lALDX)
     &       +          rk(110)*yh(lGLYD)
     &       +          rk(113)*yh(lGLY)
     &       +          rk(118)*yh(lMGLY)
     &       +          rk(120)*yh(lACET)
     &       +          rk(122)*yh(lKET)
     &       +          rk(123)*yh(lHACT)
     &       +          rk(124)*yh(lFACD)
     &       +          rk(125)*yh(lAACD)
     &       +          rk(126)*yh(lPACD)
     &       +          rk(127)*CH4
     &       +          rk(128)*yh(lECH4)
     &       +          rk(129)*yh(lETHA)
     &       +          rk(130)*yh(lPRPA)
     &       +          rk(133)*yh(lPAR)
     &       + ( 0.300)*rk(138)*yh(lETHY)
     &       +          rk(139)*yh(lETH)
     &       +          rk(142)*yh(lOLE)
     &       +          rk(145)*yh(lIOLE)
     &       + ( 0.882)*rk(148)*yh(lBENZ)
     &       + ( 0.900)*rk(152)*yh(lTOL)
     &       + ( 0.756)*rk(156)*yh(lXYL)
     &       +          rk(161)*yh(lOPEN)
     &       +          rk(165)*yh(lXOPN)
     &       +          rk(168)*yh(lCRES)
     &       +          rk(172)*yh(lCRON)
     &       +          rk(175)*yh(lCAT1)
     &       + ( 0.500)*rk(177)*yh(lARPX)
     &       +          rk(178)*yh(lISOP)
     &       + ( 0.900)*rk(185)*yh(lISPD)
     &       + (-0.100)*rk(190)*yh(lHPLD)
     &       +          rk(191)*yh(lEPOX)
     &       +          rk(195)*yh(lINTR)
     &       +          rk(196)*yh(lAPIN)
     &       +          rk(202)*yh(lTERP)
     &       +          rk(208)*yh(lSQT)
     &       +          rk(223)*yh(lOIO) )
c
c --- Loss of HO2, excluding self-reaction
c     (net either HO2 or OH produced)
c
      lsHO2 = 1.0 + rHO2_OH + rHO2_PNA + dt*(
     &       +          rk( 18)*yh(lOH)
     &       + ( 0.500)*rk( 60)*yh(lC2O3)
     &       + ( 0.500)*rk( 67)*yh(lCXO3)
     &       + ( 0.500)*rk( 73)*yh(lOPO3)
     &       +          rk( 79)*yh(lMEO2)
     &       +          rk( 85)*yh(lXO2H)
     &       +          rk( 88)*yh(lXO2)
     &       +          rk( 91)*yh(lXO2N)
     &       +          rk(150)*yh(lBZO2)
     &       +          rk(154)*yh(lTO2)
     &       +          rk(158)*yh(lXLO2)
     &       +          rk(171)*yh(lCRO)
     &       + ( 0.880)*rk(180)*yh(lISO2)
     &       + ( 0.170)*rk(198)*yh(lAPO2)
     &       + ( 0.860)*rk(204)*yh(lTPO2)
cgwilson-CB7     &       + ( 0.800)*rk(210)*yh(lSQO2)
     &       +          rk(219)*yh(lIO) )
c
c --- HO2 self-reaction
c
      self = dt*2.0*( 
     &       +          rk( 19)
     &       +          rk( 20)*H2O )
c
c --- Loss of HONO
c
      lsHONO = 1.0 + dt*(
     &       +          rk( 39)
     &       +          rk( 40)*yh(lOH) )
c
c --- Loss of PNA
c
      lsPNA = 1.0 + dt*(
     &       +          rk( 46)
     &       +          rk( 47)
     &       +          rk( 48)*yh(lOH) )
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
      y1(lHO2)  = MAX(1.0E-18, MAX(Q/A ,-C/Q) )
      y1(lOH)   = MAX(1.0E-18, ( ( newOH + rHO2_OH*y1(lHO2) )*lsHONO + 
     &                                        rHONO_OH*newHONO ) * t1 )
      y1(lPNA)  = MAX(1.0E-15, ( newPNA + rHO2_PNA*y1(lHO2) ) / lsPNA )
      y1(lHONO) = MAX(1.0E-15, ( newHONO + rOH_HONO*y1(lOH) ) / lsHONO )
c
      return
      end

