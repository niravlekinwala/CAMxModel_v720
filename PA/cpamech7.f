      subroutine cpamech7(r, rk, dtfact, nr, pa, npa, npa_init, ldark )

      use filunit
      use tracer
      use procan
      implicit none
c
c----CAMx v7.20 220430
c
c     CPAMECH7 computes chemical process analysis (CPA) parameters
c
c      Copyright 1996 - 2022
c     Ramboll
c     Created by the CMC version 5.2.6
c
c --- Subroutines Called:
c        none
c
c --- Called by:
c        PASETUP
c        CHEMDRIV
c
c --- Argument definitions:
c        r        - reaction rates (ppm/hr)
c        rk       - rate constants (ppm-n hr-1)
c        dtfact   - ratio of time step to averaging interval
c        nr       - number of reactions
c        pa       - output cpa parameter values
c        npa      - dimension of cpa parameter arrays
c        npa_init - number of cpa parameters calculated below
c        ldark    - flag set true at night
c
c --- Includes:
      include "camx.prm"
c
c --- Arguments:
c
      integer  nr, npa, npa_init
      real     r(nr), rk(nr), pa(MXCPA)
      real     dtfact
      logical  ldark
c
c --- Parameters:
c
c     Cut-point between NOx and VOC sensitive O3 production
c     Default = 0.35.  Recommended range is 0.15 to 0.35
      real     ratio_cut
      parameter (ratio_cut = 0.35)
c
c     Convert PPM to PPB
      real     ppbfact
      parameter (ppbfact = 1000.0)
c
c --- Local variables:
c
      integer  n, nn
      real     sum, sun, po3, do3, HO2toNO2, HO2_loss
      real     NO2wOH, HO2wHO2, ratio_ind
      real     OH_new, HO2_new, OH_loss, RO2_loss, HOx_CL
      real     newOH_O1D, newOH_O3, newOH_phot
      real     newHO2_O3, newHO2_pht
      real     NO3wVOC, N2O5wH2O, PAN_prdNet, ON_prod
c
c --- Entry point:
c
      nn = 0
      sun = 1.0
      if (ldark) sun = 0.0
c
c --- J(NO2) photolysis rate (per hour)
c
      nn = nn + 1
      ptname(nn)  = 'J_NO2'
      cpadesc(nn) = 'J(NO2) photolysis rate'
      cpaunit(nn) = 'hr-1'
      PA(nn) =      rk(  1)*dtfact
c
c --- J(O3O1D) photolysis rate (per hour)
c
      nn = nn + 1
      ptname(nn)  = 'J_O3O1D'
      cpadesc(nn) = 'J(O3) to O(1D) photolysis rate'
      cpaunit(nn) = 'hr-1'
      PA(nn) =      rk(  9)*dtfact
c
c
c --- Net O3 production
c
      nn = nn + 1
      ptname(nn)  = 'PO3_net'
      cpadesc(nn) = 'Net O3 production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  2)  ! O2 + O + M =
     &       + ( 0.130)*r( 60)  ! C2O3 + HO2 =
     &       + ( 0.060)*r( 67)  ! CXO3 + HO2 =
     &       + ( 0.130)*r( 73)  ! OPO3 + HO2 =
     &       -          r(  3)  ! NO + O3 =
     &       -          r(  7)  ! O3 + O =
     &       -          r(  8)  ! O3 =
     &       -          r(  9)  ! O3 =
     &       -          r( 12)  ! O3 + OH =
     &       -          r( 13)  ! O3 + HO2 =
     &       -          r( 26)  ! NO2 + O3 =
     &       -          r(140)  ! ETH + O3 =
     &       -          r(143)  ! OLE + O3 =
     &       -          r(146)  ! IOLE + O3 =
     &       -          r(162)  ! OPEN + O3 =
     &       -          r(166)  ! XOPN + O3 =
     &       -          r(183)  ! ISOP + O3 =
     &       -          r(200)  ! APIN + O3 =
     &       -          r(206)  ! TERP + O3 =
     &       -          r(212)  ! SQT + O3 =
     &       -          r(216)  ! I + O3 =
c
      po3 = max(0.0, PA(nn))
      do3 = min(0.0, PA(nn))
c
c --- Calculate the P(H2O2)/P(HNO3) indicator ratio and apply to PO3
c
      HO2wHO2 = + r( 19) + r( 20)
      NO2wOH = + r( 41) + r( 42)
      ratio_ind = min( 10., HO2wHO2/max( 1.0E-12, NO2wOH ) )*sun
c
      nn =  nn + 2
      ptname(nn-1)  = 'PO3_VOCsns'
      cpadesc(nn-1) = 'Net O3 production rate under VOC-limited condition'
      cpaunit(nn-1) = 'ppb hr-1'
      ptname(nn)  = 'PO3_NOxsns'
      cpadesc(nn) = 'Net O3 production rate under NOx-limited condition'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn-1) = 0.0
      PA(nn)   = 0.0
      if (ratio_ind .GT. 0.0) then
        if (ratio_ind .LT. ratio_cut ) then
          PA(nn-1) = PO3
        else
          PA(nn)   = PO3
        endif
      endif
c
c --- Report the P(H2O2)/P(HNO3) indicator ratio
c
      if(.NOT. lcpacum) then
        nn =  nn + 1
        ptname(nn)   = 'PH2O2_PHN3'
        cpadesc(nn) = 'Ratio of production rates for H2O2/HNO3'
        cpaunit(nn) = 'dimensionless'
        if (ldark) then
          PA(nn) = 0.0
        else
          PA(nn) = ratio_ind * (dtfact/ppbfact)
        endif
      endif
c
c --- Ozone destruction calculation from OSAT
c
      nn = nn + 1
      ptname(nn)  = 'O3_dest'
      cpadesc(nn) = 'O3 destruction rate (same implementation as OSAT)'
      cpaunit(nn) = 'ppb hr-1'
c
c  +  O1D + H2O
c
      PA(nn) =
     &       - r( 11)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
      PA(nn) = PA(nn)
     &       - r( 13)  ! O3 + HO2 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
      HO2toNO2 = 
     &       +          r( 25)  ! NO + HO2 =
     &       +          r( 32)  ! NO3 + HO2 =
c
      HO2_loss = 
     &       +          r( 13)  ! O3 + HO2 =
     &       +          r( 15)  ! HO2 + O =
     &       +          r( 18)  ! OH + HO2 =
     &       + ( 2.000)*r( 19)  ! HO2 + HO2 =
     &       + ( 2.000)*r( 20)  ! HO2 + HO2 + H2O =
     &       +          r( 25)  ! NO + HO2 =
     &       +          r( 32)  ! NO3 + HO2 =
     &       +          r( 60)  ! C2O3 + HO2 =
     &       + ( 0.750)*r( 67)  ! CXO3 + HO2 =
     &       +          r( 73)  ! OPO3 + HO2 =
     &       +          r( 79)  ! MEO2 + HO2 =
     &       +          r( 85)  ! XO2H + HO2 =
     &       +          r( 88)  ! XO2 + HO2 =
     &       +          r( 91)  ! XO2N + HO2 =
     &       +          r(150)  ! BZO2 + HO2 =
     &       +          r(154)  ! TO2 + HO2 =
     &       +          r(158)  ! XLO2 + HO2 =
     &       +          r(171)  ! CRO + HO2 =
     &       + ( 0.940)*r(180)  ! ISO2 + HO2 =
     &       + ( 0.520)*r(198)  ! APO2 + HO2 =
     &       + ( 0.930)*r(204)  ! TPO2 + HO2 =
     &       + ( 0.900)*r(210)  ! SQO2 + HO2 =
     &       +          r(219)  ! IO + HO2 =
c
      PA(nn) = PA(nn) - (
     &       + r( 12)  ! O3 + OH =
     &                     ) * (HO2_loss-HO2toNO2)/max( 1.0E-12, HO2_loss )
c
c  +  O3 + VOC
c
      PA(nn) = PA(nn)
     &       - r(140)  ! ETH + O3 =
     &       - r(143)  ! OLE + O3 =
     &       - r(146)  ! IOLE + O3 =
     &       - r(162)  ! OPEN + O3 =
     &       - r(166)  ! XOPN + O3 =
     &       - r(183)  ! ISOP + O3 =
     &       - r(200)  ! APIN + O3 =
     &       - r(206)  ! TERP + O3 =
     &       - r(212)  ! SQT + O3 =
c
c  +  O(3P) + VOC
c
      PA(nn) = PA(nn)
c
c  +  Halogen catalytic destruction (the limiting reactions)
c
      PA(nn) = PA(nn)
     &       -          r(219)  ! IO + HO2 =
     &       - ( 0.400)*r(218)  ! IO + IO =
     &       -          r(221)  ! IO + NO2 =
c
      PA(nn) = min(PA(nn), do3)
c
c --- Total OH production
c
      nn = nn + 1
      ptname(nn)  = 'OH_prod'
      cpadesc(nn) = 'Total OH production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 11)  ! O1D + H2O =
     &       +          r( 13)  ! O3 + HO2 =
     &       +          r( 15)  ! HO2 + O =
     &       + ( 2.000)*r( 21)  ! H2O2 =
     &       +          r( 23)  ! H2O2 + O =
     &       +          r( 25)  ! NO + HO2 =
     &       +          r( 32)  ! NO3 + HO2 =
     &       +          r( 39)  ! HONO =
     &       +          r( 44)  ! HNO3 =
     &       + ( 0.410)*r( 47)  ! PNA =
     &       + ( 0.500)*r( 60)  ! C2O3 + HO2 =
     &       + ( 0.250)*r( 67)  ! CXO3 + HO2 =
     &       + ( 0.500)*r( 73)  ! OPO3 + HO2 =
     &       +          r( 83)  ! MEPX =
     &       +          r( 94)  ! ROOH =
     &       + ( 0.190)*r(111)  ! GLYD =
     &       + ( 0.170)*r(140)  ! ETH + O3 =
     &       + ( 0.334)*r(143)  ! OLE + O3 =
     &       + ( 0.500)*r(146)  ! IOLE + O3 =
     &       + ( 0.500)*r(162)  ! OPEN + O3 =
     &       + ( 0.500)*r(166)  ! XOPN + O3 =
     &       + ( 0.060)*r(180)  ! ISO2 + HO2 =
     &       + ( 1.700)*r(182)  ! ISO2 =
     &       + ( 0.280)*r(183)  ! ISOP + O3 =
     &       + ( 0.250)*r(184)  ! ISOP + NO3 =
     &       + ( 0.200)*r(187)  ! ISPD =
     &       + ( 1.650)*r(189)  ! HPLD =
     &       + ( 0.100)*r(190)  ! HPLD + OH =
     &       + ( 0.700)*r(192)  ! EPX2 + NO =
     &       + ( 1.700)*r(193)  ! EPX2 + HO2 =
     &       + ( 0.850)*r(194)  ! EPX2 + RO2 =
     &       + ( 0.350)*r(198)  ! APO2 + HO2 =
     &       + ( 0.770)*r(200)  ! APIN + O3 =
     &       + ( 0.700)*r(201)  ! APIN + NO3 =
     &       + ( 0.070)*r(204)  ! TPO2 + HO2 =
     &       + ( 0.470)*r(206)  ! TERP + O3 =
     &       + ( 0.480)*r(207)  ! TERP + NO3 =
     &       + ( 0.100)*r(210)  ! SQO2 + HO2 =
     &       + ( 0.080)*r(212)  ! SQT + O3 =
     &       + ( 0.470)*r(213)  ! SQT + NO3 =
     &       +          r(215)  ! HOI =
c
c --- OH from O(1D)
c
      nn = nn + 1
      ptname(nn)  = 'newOH_O1D'
      cpadesc(nn) = 'OH production rate from O(1D) + H2O'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = ( 2.000)*r( 11)  ! O1D + H2O =
      newOH_O1D = PA(nn)
c
c --- OH from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newOH_O3'
      cpadesc(nn) = 'OH production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.170)*r(140)  ! ETH + O3 =
     &       + ( 0.334)*r(143)  ! OLE + O3 =
     &       + ( 0.500)*r(146)  ! IOLE + O3 =
     &       + ( 0.500)*r(162)  ! OPEN + O3 =
     &       + ( 0.500)*r(166)  ! XOPN + O3 =
     &       + ( 0.280)*r(183)  ! ISOP + O3 =
     &       + ( 0.770)*r(200)  ! APIN + O3 =
     &       + ( 0.470)*r(206)  ! TERP + O3 =
     &       + ( 0.080)*r(212)  ! SQT + O3 =
      newOH_O3 = PA(nn)
c
c --- OH directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newOH_phot'
      cpadesc(nn) = 'OH production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 21)  ! H2O2 =
     &       +          r( 39)  ! HONO =
     &       +          r( 44)  ! HNO3 =
     &       +          r( 83)  ! MEPX =
     &       +          r( 94)  ! ROOH =
     &       + ( 0.190)*r(111)  ! GLYD =
     &       + ( 1.700)*r(182)  ! ISO2 =
     &       + ( 0.200)*r(187)  ! ISPD =
     &       + ( 1.650)*r(189)  ! HPLD =
      newOH_phot = PA(nn)
c
c --- New OH
c
      nn = nn + 1
      ptname(nn)  = 'OH_new'
      cpadesc(nn) = 'New OH, i.e., newOH_O1D + newOH_O3 + newOH_phot'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = newOH_O1D + newOH_O3 + newOH_phot
      OH_new = PA(nn)
c
c --- OH from HONO
c
      nn = nn + 1
      ptname(nn)  = 'newOH_HONO'
      cpadesc(nn) = 'OH production rate from HONO photolysis'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =          r( 39)  ! HONO =
c
c --- OH from HPLD
c
      nn = nn + 1
      ptname(nn)  = 'newOH_HPLD'
      cpadesc(nn) = 'OH production rate from HPALD photolysis'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = ( 1.650)*r(189)  ! HPLD = ( 0.100)*r(190)  ! HPLD + OH = ( 0.100)*r(190)  ! HPLD + OH =
c
c --- OH loss
c
      nn = nn + 1
      ptname(nn)  = 'OH_loss'
      cpadesc(nn) = 'OH loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 12)  ! O3 + OH =
     &       +          r( 14)  ! OH + O =
     &       + ( 2.000)*r( 16)  ! OH + OH =
     &       + ( 2.000)*r( 17)  ! OH + OH =
     &       +          r( 18)  ! OH + HO2 =
     &       +          r( 22)  ! H2O2 + OH =
     &       +          r( 31)  ! NO3 + OH =
     &       +          r( 38)  ! NO + OH =
     &       +          r( 40)  ! HONO + OH =
     &       +          r( 41)  ! NO2 + OH =
     &       +          r( 42)  ! NO2 + OH + H2O =
     &       +          r( 43)  ! HNO3 + OH =
     &       +          r( 48)  ! PNA + OH =
     &       +          r( 49)  ! H2 + OH =
     &       +          r( 50)  ! CO + OH =
     &       +          r( 51)  ! SO2 + OH =
     &       +          r( 53)  ! DMS + OH =
     &       +          r( 54)  ! DMS + OH + O2 =
     &       +          r( 66)  ! PANX + OH =
     &       +          r( 72)  ! OPAN + OH =
     &       + ( 0.600)*r( 82)  ! MEPX + OH =
     &       + ( 0.600)*r( 93)  ! ROOH + OH =
     &       +          r( 95)  ! NTR1 + OH =
     &       +          r( 98)  ! MEOH + OH =
     &       +          r( 99)  ! ETOH + OH =
     &       +          r(100)  ! FORM + OH =
     &       +          r(104)  ! ALD2 + OH =
     &       +          r(107)  ! ALDX + OH =
     &       +          r(110)  ! GLYD + OH =
     &       +          r(113)  ! GLY + OH =
     &       +          r(118)  ! MGLY + OH =
     &       +          r(120)  ! ACET + OH =
     &       +          r(122)  ! KET + OH =
     &       +          r(123)  ! HACT + OH =
     &       +          r(124)  ! FACD + OH =
     &       +          r(125)  ! AACD + OH =
     &       +          r(126)  ! PACD + OH =
     &       +          r(127)  ! CH4 + OH =
     &       +          r(128)  ! ECH4 + OH =
     &       +          r(129)  ! ETHA + OH =
     &       +          r(130)  ! PRPA + OH =
     &       +          r(133)  ! PAR + OH =
     &       + ( 0.300)*r(138)  ! ETHY + OH =
     &       +          r(139)  ! ETH + OH =
     &       +          r(142)  ! OLE + OH =
     &       +          r(145)  ! IOLE + OH =
     &       + ( 0.882)*r(148)  ! BENZ + OH =
     &       + ( 0.900)*r(152)  ! TOL + OH =
     &       + ( 0.756)*r(156)  ! XYL + OH =
     &       +          r(161)  ! OPEN + OH =
     &       +          r(165)  ! XOPN + OH =
     &       +          r(168)  ! CRES + OH =
     &       +          r(172)  ! CRON + OH =
     &       +          r(175)  ! CAT1 + OH =
     &       + ( 0.500)*r(177)  ! ARPX + OH =
     &       +          r(178)  ! ISOP + OH =
     &       + ( 0.900)*r(185)  ! ISPD + OH =
     &       +          r(191)  ! EPOX + OH =
     &       +          r(195)  ! INTR + OH =
     &       +          r(196)  ! APIN + OH =
     &       +          r(202)  ! TERP + OH =
     &       +          r(208)  ! SQT + OH =
     &       +          r(223)  ! OIO + OH =
      OH_loss = PA(nn)
c
c --- OH with CO
c
      nn = nn + 1
      ptname(nn)  = 'OHwCO'
      cpadesc(nn) = 'OH reaction rate with CO'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 50)  ! CO + OH =
c
c --- OH with ECH4
c
      nn = nn + 1
      ptname(nn)  = 'OHwECH4'
      cpadesc(nn) = 'OH reaction rate with ECH4 (emitted CH4)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(128)  ! ECH4 + OH =
c
c --- OH with isoprene
c
      nn = nn + 1
      ptname(nn)  = 'OHwISOP'
      cpadesc(nn) = 'OH reaction rate with isoprene'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(178)  ! ISOP + OH =
c
c --- OH with VOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwVOC'
      cpadesc(nn) = 'OH reaction rate with all VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 98)  ! MEOH + OH =
     &       +          r( 99)  ! ETOH + OH =
     &       +          r(100)  ! FORM + OH =
     &       +          r(104)  ! ALD2 + OH =
     &       +          r(107)  ! ALDX + OH =
     &       +          r(120)  ! ACET + OH =
     &       +          r(122)  ! KET + OH =
     &       +          r(129)  ! ETHA + OH =
     &       +          r(130)  ! PRPA + OH =
     &       +          r(133)  ! PAR + OH =
     &       +          r(138)  ! ETHY + OH =
     &       +          r(139)  ! ETH + OH =
     &       +          r(142)  ! OLE + OH =
     &       +          r(145)  ! IOLE + OH =
     &       +          r(148)  ! BENZ + OH =
     &       +          r(152)  ! TOL + OH =
     &       +          r(156)  ! XYL + OH =
     &       +          r(178)  ! ISOP + OH =
     &       +          r(196)  ! APIN + OH =
     &       +          r(202)  ! TERP + OH =
     &       +          r(208)  ! SQT + OH =
c
c --- OH with HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwHRVOC'
      cpadesc(nn) = 'OH reaction rate with HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(139)  ! ETH + OH =
     &       +          r(142)  ! OLE + OH =
     &       +          r(145)  ! IOLE + OH =
c
c --- OH with Aromatics
c
      nn = nn + 1
      ptname(nn)  = 'OHwArom'
      cpadesc(nn) = 'OH reaction rate with aromatic VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(148)  ! BENZ + OH =
     &       +          r(152)  ! TOL + OH =
     &       +          r(156)  ! XYL + OH =
c
c --- OH with Alkanes (except methane)
c
      nn = nn + 1
      ptname(nn)  = 'OHwAlkane'
      cpadesc(nn) = 'OH reaction rate with alkanes (except methane)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(129)  ! ETHA + OH =
     &       +          r(130)  ! PRPA + OH =
     &       +          r(133)  ! PAR + OH =
c
c --- Total HCHO production
c
C      nn = nn + 1
C      ptname(nn)  = 'HCHO_prod'
C      cpadesc(nn) = 'Total HCHO production rate'
C      cpaunit(nn) = 'ppb hr-1'
C      PA(nn) =
C     &       + ( 0.200)*r(190)  ! HPLD + OH =
c
c --- HCHO production from HRVOC
c
C      nn = nn + 1
C      ptname(nn)  = 'nwHCHO_HRV'
C      cpadesc(nn) = 'HCHO production rate from HRVOC'
C      cpaunit(nn) = 'ppb hr-1'
C      PA(nn) =
c
c --- HCHO production from Isoprene at first generation
c
C      nn = nn + 1
C      ptname(nn)  = 'nwHCHO_ISP'
C      cpadesc(nn) = 'HCHO production rate from isoprene at first generation'
C      cpaunit(nn) = 'ppb hr-1'
C      PA(nn) =
c
c --- Total HO2 production, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'HO2_prod'
      cpadesc(nn) = 'Total HO2 production rate, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 12)  ! O3 + OH =
     &       +          r( 14)  ! OH + O =
     &       +          r( 22)  ! H2O2 + OH =
     &       +          r( 23)  ! H2O2 + O =
     &       +          r( 31)  ! NO3 + OH =
     &       +          r( 49)  ! H2 + OH =
     &       +          r( 50)  ! CO + OH =
     &       +          r( 51)  ! SO2 + OH =
     &       + ( 0.350)*r( 68)  ! CXO3 + RO2 =
     &       + ( 0.800)*r( 69)  ! OPO3 + NO =
     &       +          r( 78)  ! MEO2 + NO =
     &       + ( 0.900)*r( 80)  ! MEO2 + C2O3 =
     &       + ( 0.370)*r( 81)  ! MEO2 + RO2 =
     &       +          r( 84)  ! XO2H + NO =
     &       + ( 0.600)*r( 86)  ! XO2H + RO2 =
     &       +          r( 94)  ! ROOH =
     &       +          r( 98)  ! MEOH + OH =
     &       + ( 0.900)*r( 99)  ! ETOH + OH =
     &       +          r(100)  ! FORM + OH =
     &       + ( 2.000)*r(101)  ! FORM =
     &       +          r(103)  ! FORM + NO3 =
     &       +          r(106)  ! ALD2 =
     &       +          r(109)  ! ALDX =
     &       + ( 0.200)*r(110)  ! GLYD + OH =
     &       + ( 1.400)*r(111)  ! GLYD =
     &       +          r(113)  ! GLY + OH =
     &       + ( 2.000)*r(114)  ! GLY =
     &       +          r(115)  ! GLY + NO3 =
     &       +          r(116)  ! MGLY =
     &       + ( 0.960)*r(121)  ! KET =
     &       + ( 0.610)*r(122)  ! KET + OH =
     &       +          r(123)  ! HACT + OH =
     &       +          r(124)  ! FACD + OH =
     &       + ( 0.620)*r(136)  ! ROR =
     &       +          r(137)  ! ROR + O2 =
     &       + ( 0.300)*r(138)  ! ETHY + OH =
     &       + ( 0.270)*r(140)  ! ETH + O3 =
     &       + ( 0.080)*r(143)  ! OLE + O3 =
     &       + ( 0.530)*r(148)  ! BENZ + OH =
     &       + ( 0.918)*r(149)  ! BZO2 + NO =
     &       +          r(151)  ! BZO2 + RO2 =
     &       + ( 0.180)*r(152)  ! TOL + OH =
     &       + ( 0.860)*r(153)  ! TO2 + NO =
     &       +          r(155)  ! TO2 + RO2 =
     &       + ( 0.155)*r(156)  ! XYL + OH =
     &       + ( 0.860)*r(157)  ! XLO2 + NO =
     &       +          r(159)  ! XLO2 + RO2 =
     &       +          r(160)  ! OPEN =
     &       + ( 0.560)*r(162)  ! OPEN + O3 =
     &       + ( 0.700)*r(164)  ! XOPN =
     &       +          r(168)  ! CRES + OH =
     &       +          r(174)  ! CRON =
     &       + ( 0.200)*r(175)  ! CAT1 + OH =
     &       + ( 0.900)*r(179)  ! ISO2 + NO =
     &       + ( 0.350)*r(182)  ! ISO2 =
     &       + ( 0.500)*r(183)  ! ISOP + O3 =
     &       + ( 0.100)*r(185)  ! ISPD + OH =
     &       + ( 0.500)*r(188)  ! ISPX + OH =
     &       + ( 0.200)*r(189)  ! HPLD =
     &       + ( 0.200)*r(191)  ! EPOX + OH =
     &       +          r(192)  ! EPX2 + NO =
     &       +          r(194)  ! EPX2 + RO2 =
     &       + ( 0.770)*r(197)  ! APO2 + NO =
     &       + ( 0.500)*r(199)  ! APO2 + RO2 =
     &       + ( 0.170)*r(200)  ! APIN + O3 =
     &       + ( 0.700)*r(201)  ! APIN + NO3 =
     &       + ( 0.750)*r(203)  ! TPO2 + NO =
     &       + ( 0.500)*r(205)  ! TPO2 + RO2 =
     &       + ( 0.160)*r(206)  ! TERP + O3 =
     &       + ( 0.300)*r(207)  ! TERP + NO3 =
     &       + ( 0.600)*r(209)  ! SQO2 + NO =
     &       + ( 0.500)*r(211)  ! SQO2 + RO2 =
     &       + ( 0.080)*r(212)  ! SQT + O3 =
c
c --- HO2 from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_O3'
      cpadesc(nn) = 'HO2 production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.270)*r(140)  ! ETH + O3 =
     &       + ( 0.080)*r(143)  ! OLE + O3 =
     &       + ( 0.560)*r(162)  ! OPEN + O3 =
     &       + ( 0.500)*r(183)  ! ISOP + O3 =
     &       + ( 0.170)*r(200)  ! APIN + O3 =
     &       + ( 0.160)*r(206)  ! TERP + O3 =
     &       + ( 0.080)*r(212)  ! SQT + O3 =
      newHO2_O3 = PA(nn)
c
c --- HO2 directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_pht'
      cpadesc(nn) = 'HO2 production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 94)  ! ROOH =
     &       + ( 2.000)*r(101)  ! FORM =
     &       +          r(106)  ! ALD2 =
     &       +          r(109)  ! ALDX =
     &       + ( 1.400)*r(111)  ! GLYD =
     &       + ( 2.000)*r(114)  ! GLY =
     &       +          r(116)  ! MGLY =
     &       + ( 0.960)*r(121)  ! KET =
     &       + ( 0.620)*r(136)  ! ROR =
     &       +          r(160)  ! OPEN =
     &       + ( 0.700)*r(164)  ! XOPN =
     &       +          r(174)  ! CRON =
     &       + ( 0.350)*r(182)  ! ISO2 =
     &       + ( 0.200)*r(189)  ! HPLD =
      newHO2_pht = PA(nn)
c
c --- New HO2
c
      nn = nn + 1
      ptname(nn)  = 'HO2_new'
      cpadesc(nn) = 'New HO2 production rate, i.e. newHO2_O3 + newHO2_pht'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = newHO2_O3 + newHO2_pht
      OH_new = PA(nn)
c
c --- HO2 from HCHO photolysis
c
C      nn = nn + 1
C      ptname(nn)  = 'nwHO2_HCHO'
C      cpadesc(nn) = 'HO2 production rate from HCHO photolysis'
C      cpaunit(nn) = 'ppb hr-1'
C      PA(nn) =
c
c --- HO2 loss (except to HNO4) - computed above
c
      nn = nn + 1
      ptname(nn)  = 'HO2_loss'
      cpadesc(nn) = 'HO2 loss (except to HNO4)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = HO2_loss
c
c --- HO2wHO2 - computed above
c
      nn = nn + 1
      ptname(nn)  = 'HO2wHO2'
      cpadesc(nn) = 'HO2 + HO2 self-reaction'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = HO2wHO2
c
c --- HOx chain length = HOx reacted / new HOx
c
      if(.NOT. lcpacum) then
        nn = nn + 1
        ptname(nn)  = 'HOx_CL'
        cpadesc(nn) = 'HOx chain length = HOx reacted / new HOx'
        cpaunit(nn) = 'dimensionless'
        PA(nn) = ((OH_loss + HO2_loss) / max( 1.0E-12, (OH_new + HO2_new)))
     &           *sun*(dtfact/ppbfact)
      endif
c
c --- NO3 production
c
      nn = nn + 1
      ptname(nn)  = 'NO3_prod'
      cpadesc(nn) = 'NO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  6)  ! NO2 + O =
     &       +          r( 26)  ! NO2 + O3 =
     &       +          r( 35)  ! N2O5 =
     &       +          r( 36)  ! N2O5 =
     &       +          r( 43)  ! HNO3 + OH =
     &       + ( 0.410)*r( 47)  ! PNA =
     &       + ( 0.400)*r( 59)  ! PAN =
     &       +          r(228)  ! INO3 =
c
c --- NO3 production from N2O5
c
      nn = nn + 1
      ptname(nn)  = 'N2O5toNO3'
      cpadesc(nn) = 'NO3 production from N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 35)  ! N2O5 =
     &       +          r( 36)  ! N2O5 =
c
c --- NO3 loss
c
      nn = nn + 1
      ptname(nn)  = 'NO3_loss'
      cpadesc(nn) = 'NO3 loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 27)  ! NO3 =
     &       +          r( 28)  ! NO3 =
     &       +          r( 29)  ! NO3 + NO =
     &       +          r( 30)  ! NO3 + NO2 =
     &       +          r( 31)  ! NO3 + OH =
     &       +          r( 32)  ! NO3 + HO2 =
     &       + ( 2.000)*r( 33)  ! NO3 + NO3 =
     &       +          r( 34)  ! NO3 + NO2 =
     &       +          r( 55)  ! DMS + NO3 =
     &       +          r(103)  ! FORM + NO3 =
     &       +          r(105)  ! ALD2 + NO3 =
     &       +          r(108)  ! ALDX + NO3 =
     &       +          r(112)  ! GLYD + NO3 =
     &       +          r(115)  ! GLY + NO3 =
     &       +          r(117)  ! MGLY + NO3 =
     &       +          r(141)  ! ETH + NO3 =
     &       +          r(144)  ! OLE + NO3 =
     &       +          r(147)  ! IOLE + NO3 =
     &       +          r(163)  ! OPEN + NO3 =
     &       +          r(167)  ! XOPN + NO3 =
     &       +          r(169)  ! CRES + NO3 =
     &       +          r(173)  ! CRON + NO3 =
     &       +          r(176)  ! CAT1 + NO3 =
     &       +          r(184)  ! ISOP + NO3 =
     &       +          r(186)  ! ISPD + NO3 =
     &       +          r(201)  ! APIN + NO3 =
     &       +          r(207)  ! TERP + NO3 =
     &       +          r(213)  ! SQT + NO3 =
c
c --- NO3 to N2O5
c
      nn = nn + 1
      ptname(nn)  = 'NO3toN2O5'
      cpadesc(nn) = 'NO3 conversion to N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 34)  ! NO3 + NO2 =
c
c --- RO2 loss
c
      nn = nn + 1
      ptname(nn)  = 'RO2_loss'
      cpadesc(nn) = 'RO2 loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 75)  ! RO2 + NO =
     &       +          r( 76)  ! RO2 + HO2 =
     &       + ( 2.000)*r( 77)  ! RO2 + RO2 =
      RO2_loss = PA(nn)
c
c --- RO2 with NO
c
      nn = nn + 1
      ptname(nn)  = 'RO2wNO'
      cpadesc(nn) = 'RO2 reaction with NO'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 75)  ! RO2 + NO =
c
c --- RO2 with HO2
c
      nn = nn + 1
      ptname(nn)  = 'RO2wHO2'
      cpadesc(nn) = 'RO2 reaction with HO2'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 76)  ! RO2 + HO2 =
c
c --- RO2 with RO2
c
      nn = nn + 1
      ptname(nn)  = 'RO2wRO2'
      cpadesc(nn) = 'RO2 + RO2 self-reaction'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 77)  ! RO2 + RO2 =
c
c --- Total production of ON compounds
c
      nn = nn + 1
      ptname(nn)  = 'ON_prod'
      cpadesc(nn) = 'Total production of organic nitrate (ON)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.500)*r( 72)  ! OPAN + OH =
     &       + ( 0.500)*r( 90)  ! XO2N + NO =
     &       + ( 0.082)*r(149)  ! BZO2 + NO =
     &       + ( 0.140)*r(153)  ! TO2 + NO =
     &       + ( 0.140)*r(157)  ! XLO2 + NO =
     &       + ( 0.500)*r(167)  ! XOPN + NO3 =
     &       + ( 0.750)*r(184)  ! ISOP + NO3 =
     &       + ( 0.900)*r(186)  ! ISPD + NO3 =
     &       + ( 0.020)*r(192)  ! EPX2 + NO =
     &       + ( 0.230)*r(197)  ! APO2 + NO =
     &       + ( 0.300)*r(201)  ! APIN + NO3 =
     &       + ( 0.260)*r(203)  ! TPO2 + NO =
     &       + ( 0.700)*r(207)  ! TERP + NO3 =
     &       + ( 0.400)*r(209)  ! SQO2 + NO =
     &       + ( 0.760)*r(213)  ! SQT + NO3 =
     &       + ( 0.500)*r( 90)  ! XO2N + NO =
     &       + ( 0.500)*r(141)  ! ETH + NO3 =
     &       + ( 0.500)*r(144)  ! OLE + NO3 =
     &       + ( 0.500)*r(147)  ! IOLE + NO3 =
     &       +          r(170)  ! CRO + NO2 =
     &       + ( 0.100)*r(179)  ! ISO2 + NO =
      ON_prod = PA(nn)
c
c --- INTR production
c
      nn = nn + 1
      ptname(nn)  = 'INTR_prod'
      cpadesc(nn) = 'Total production of the species INTR'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.100)*r(179)  ! ISO2 + NO =
c
c --- NTR1 production
c
      nn = nn + 1
      ptname(nn)  = 'NTR1_prod'
      cpadesc(nn) = 'Total production of the species NTR1'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.500)*r( 90)  ! XO2N + NO =
     &       + ( 0.500)*r(141)  ! ETH + NO3 =
     &       + ( 0.500)*r(144)  ! OLE + NO3 =
     &       + ( 0.500)*r(147)  ! IOLE + NO3 =
c
c --- NTR2 production
c
      nn = nn + 1
      ptname(nn)  = 'NTR2_prod'
      cpadesc(nn) = 'Total production of the species NTR2'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.500)*r( 72)  ! OPAN + OH =
     &       + ( 0.500)*r( 90)  ! XO2N + NO =
     &       + ( 0.082)*r(149)  ! BZO2 + NO =
     &       + ( 0.140)*r(153)  ! TO2 + NO =
     &       + ( 0.140)*r(157)  ! XLO2 + NO =
     &       + ( 0.500)*r(167)  ! XOPN + NO3 =
     &       +          r(172)  ! CRON + OH =
     &       +          r(173)  ! CRON + NO3 =
     &       + ( 0.750)*r(184)  ! ISOP + NO3 =
     &       + ( 0.900)*r(186)  ! ISPD + NO3 =
     &       + ( 0.020)*r(192)  ! EPX2 + NO =
     &       + ( 0.400)*r(195)  ! INTR + OH =
     &       + ( 0.230)*r(197)  ! APO2 + NO =
     &       + ( 0.300)*r(201)  ! APIN + NO3 =
     &       + ( 0.260)*r(203)  ! TPO2 + NO =
     &       + ( 0.700)*r(207)  ! TERP + NO3 =
     &       + ( 0.400)*r(209)  ! SQO2 + NO =
     &       + ( 0.760)*r(213)  ! SQT + NO3 =
c
c --- NTR2 production from NTR1 with OH
c
C      nn = nn + 1
C      ptname(nn)  = 'NTR1wOH'
C      cpadesc(nn) = 'NTR2 production from NTR1 + OH'
C      cpaunit(nn) = 'ppb hr-1'
C      PA(nn) =
c
c --- HNO3 production
c
      nn = nn + 1
      ptname(nn)  = 'HNO3_prod'
      cpadesc(nn) = 'HNO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 37)  ! N2O5 + H2O =
     &       +          r( 41)  ! NO2 + OH =
     &       +          r( 42)  ! NO2 + OH + H2O =
     &       +          r( 55)  ! DMS + NO3 =
     &       +          r( 97)  ! NTR2 =
     &       +          r(103)  ! FORM + NO3 =
     &       +          r(105)  ! ALD2 + NO3 =
     &       +          r(108)  ! ALDX + NO3 =
     &       +          r(112)  ! GLYD + NO3 =
     &       +          r(115)  ! GLY + NO3 =
     &       +          r(117)  ! MGLY + NO3 =
     &       +          r(163)  ! OPEN + NO3 =
     &       +          r(169)  ! CRES + NO3 =
     &       +          r(173)  ! CRON + NO3 =
     &       +          r(176)  ! CAT1 + NO3 =
     &       + ( 0.100)*r(186)  ! ISPD + NO3 =
     &       +          r(229)  ! INO3 + H2O =
c
c --- HNO3 production from NO2 with OH - computed above
c
      nn = nn + 1
      ptname(nn)  = 'NO2wOH'
      cpadesc(nn) = 'HNO3 production from OH + NO2'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = NO2wOH
c
c --- HNO3 production from N2O5 with H2O
c
      nn = nn + 1
      ptname(nn)  = 'N2O5wH2O'
      cpadesc(nn) = 'HNO3 production from N2O5 + H2O'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 37)  ! N2O5 + H2O =
      N2O5wH2O = PA(nn)
c
c --- NO2 recycled from HNO3 and organic nitrates
c
c     explicit isoprene nitrates and nitrophenols excluded
c
      nn = nn + 1
      ptname(nn)  = 'NO2_rcycl'
      cpadesc(nn) = 'NO2 recycled from HNO3 and organic nitrates'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 44)  ! HNO3 =
     &       +          r( 95)  ! NTR1 + OH =
     &       +          r( 96)  ! NTR1 =
c
c --- Net production of PAN compounds
c
      nn = nn + 1
      ptname(nn)  = 'PAN_prdNet'
      cpadesc(nn) = 'Net production of PAN species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 57)  ! C2O3 + NO2 =
     &       +          r( 64)  ! CXO3 + NO2 =
     &       +          r( 70)  ! OPO3 + NO2 =
     &       -          r( 58)  ! PAN =
     &       -          r( 59)  ! PAN =
     &       -          r( 65)  ! PANX =
     &       -          r( 66)  ! PANX + OH =
     &       -          r( 71)  ! OPAN =
     &       -          r( 72)  ! OPAN + OH =
      PAN_prdNet = PA(nn)
c
c --- NO3 with VOC
c
      nn = nn + 1
      ptname(nn)  = 'NO3wVOC'
      cpadesc(nn) = 'NO3 reaction with VOC species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(103)  ! FORM + NO3 =
     &       +          r(105)  ! ALD2 + NO3 =
     &       +          r(108)  ! ALDX + NO3 =
     &       +          r(141)  ! ETH + NO3 =
     &       +          r(144)  ! OLE + NO3 =
     &       +          r(147)  ! IOLE + NO3 =
     &       +          r(184)  ! ISOP + NO3 =
     &       +          r(201)  ! APIN + NO3 =
     &       +          r(207)  ! TERP + NO3 =
     &       +          r(213)  ! SQT + NO3 =
      NO3wVOC = PA(nn)
c
c --- LN/Q calculation
c
      if(.NOT. lcpacum) then
        nn = nn + 1
        ptname(nn)  = 'LNoQ'
        cpadesc(nn) = 'LN/Q metric for ozone sensitivity'
        cpaunit(nn) = 'dimensionless'
        PA(nn) = ( NO2wOH + NO3wVOC + N2O5wH2O +
     &             max( 0., PAN_prdNet ) + ON_prod )
     &           / max( 1.0E-12, (OH_new + HO2_new) )
        PA(nn) = max( 0., min( 10., PA(nn) ))*sun*(dtfact/ppbfact)
      endif
c
c --- I production
c
      nn = nn + 1
      ptname(nn)  = 'I_prod'
      cpadesc(nn) = 'I atom production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r(214)  ! I2 =
     &       +          r(215)  ! HOI =
     &       +          r(217)  ! IO =
     &       + ( 0.400)*r(218)  ! IO + IO =
     &       +          r(220)  ! IO + NO =
     &       +          r(222)  ! OIO =
     &       +          r(226)  ! I2O2 =
     &       +          r(228)  ! INO3 =
c
c --- Ozone destruction by I
c
      nn = nn + 1
      ptname(nn)  = 'I_O3dest'
      cpadesc(nn) = 'Ozone destruction by I atom'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       -          r(219)  ! IO + HO2 =
     &       - ( 0.400)*r(218)  ! IO + IO =
     &       -          r(221)  ! IO + NO2 =
c
c --- Convert to ppb/hr
c
      Do n = 1,nn
        PA(n) = ppbfact*PA(n)
      End do
c
c --- Set num of outputs on first dummy call from pasetup
c
      npa_init = nn
c
      return
      end

