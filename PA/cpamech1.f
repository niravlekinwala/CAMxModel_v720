      subroutine cpamech1(r, rk, dtfact, nr, pa, npa, npa_init, ldark )

      use filunit
      use tracer
      use procan
      implicit none
c
c----CAMx v7.20 220430
c
c     CPAMECH1 computes chemical process analysis (CPA) parameters
c
c      Copyright 1996 - 2022
c     Ramboll
c     Created by the CMC version 5.2.6 on 2020 02 20
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
     &       +          r(  2)  ! O + O2 + M =
     &       + ( 0.130)*r( 57)  ! C2O3 + HO2 =
     &       + ( 0.130)*r( 65)  ! CXO3 + HO2 =
     &       + ( 0.130)*r(202)  ! OPO3 + HO2 =
     &       -          r(  3)  ! O3 + NO =
     &       -          r(  7)  ! O + O3 =
     &       -          r(  8)  ! O3 =
     &       -          r(  9)  ! O3 =
     &       -          r( 12)  ! O3 + OH =
     &       -          r( 13)  ! O3 + HO2 =
     &       -          r( 26)  ! NO2 + O3 =
     &       -          r( 34)  ! NO3 + O3 =
     &       -          r(135)  ! ETH + O3 =
     &       -          r(138)  ! OLE + O3 =
     &       -          r(141)  ! IOLE + O3 =
     &       -          r(149)  ! ISOP + O3 =
     &       -          r(152)  ! ISPD + O3 =
     &       -          r(165)  ! TERP + O3 =
     &       -          r(191)  ! XOPN + O3 =
     &       -          r(195)  ! OPEN + O3 =
     &       -          r(211)  ! I + O3 =
     &       -          r(222)  ! I2O2 + O3 =
c
      po3 = max(0.0, PA(nn))
      do3 = min(0.0, PA(nn))
c
c --- Calculate the P(H2O2)/P(HNO3) indicator ratio and apply to PO3
c
      HO2wHO2 = + r( 19) + r( 20)
      NO2wOH = + r( 45) + r(234)
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
     &       +          r( 25)  ! HO2 + NO =
     &       +          r( 33)  ! NO3 + HO2 =
c
      HO2_loss = 
     &       +          r( 13)  ! O3 + HO2 =
     &       +          r( 15)  ! HO2 + O =
     &       +          r( 18)  ! OH + HO2 =
     &       + ( 2.000)*r( 19)  ! HO2 + HO2 =
     &       + ( 2.000)*r( 20)  ! HO2 + HO2 + H2O =
     &       +          r( 25)  ! HO2 + NO =
     &       +          r( 33)  ! NO3 + HO2 =
     &       +          r( 57)  ! C2O3 + HO2 =
     &       +          r( 65)  ! CXO3 + HO2 =
     &       +          r( 72)  ! MEO2 + HO2 =
     &       +          r( 76)  ! XO2H + HO2 =
     &       +          r( 80)  ! XO2 + HO2 =
     &       +          r( 84)  ! XO2N + HO2 =
     &       +          r(100)  ! FORM + HO2 =
     &       + ( 0.800)*r(103)  ! HCO3 + HO2 =
     &       + ( 0.880)*r(145)  ! ISO2 + HO2 =
     &       + ( 0.175)*r(159)  ! EPX2 + HO2 =
     &       +          r(170)  ! BZO2 + HO2 =
     &       +          r(175)  ! TO2 + HO2 =
     &       +          r(179)  ! XLO2 + HO2 =
     &       +          r(185)  ! CRO + HO2 =
     &       +          r(202)  ! OPO3 + HO2 =
     &       +          r(214)  ! IO + HO2 =
c
      PA(nn) = PA(nn) - (
     &       + r( 12)  ! O3 + OH =
     &                     ) * (HO2_loss-HO2toNO2)/max( 1.0E-12, HO2_loss )
c
c  +  O3 + VOC
c
      PA(nn) = PA(nn)
     &       - r(135)  ! ETH + O3 =
     &       - r(138)  ! OLE + O3 =
     &       - r(141)  ! IOLE + O3 =
     &       - r(149)  ! ISOP + O3 =
     &       - r(152)  ! ISPD + O3 =
     &       - r(165)  ! TERP + O3 =
c
c  +  O(3P) + VOC
c
      PA(nn) = PA(nn)
c
c  +  Halogen catalytic destruction (the limiting reactions)
c
      PA(nn) = PA(nn)
     &       -          r(214)  ! IO + HO2 =
     &       - ( 0.400)*r(213)  ! IO + IO =
     &       -          r(216)  ! IO + NO2 =
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
     &       +          r( 25)  ! HO2 + NO =
     &       +          r( 33)  ! NO3 + HO2 =
     &       +          r( 43)  ! HONO =
     &       +          r( 47)  ! HNO3 =
     &       + ( 0.410)*r( 50)  ! PNA =
     &       + ( 0.500)*r( 57)  ! C2O3 + HO2 =
     &       + ( 0.500)*r( 65)  ! CXO3 + HO2 =
     &       +          r( 88)  ! MEPX =
     &       +          r( 90)  ! ROOH =
     &       + ( 0.200)*r(103)  ! HCO3 + HO2 =
     &       + ( 0.190)*r(111)  ! GLYD =
     &       + ( 0.170)*r(135)  ! ETH + O3 =
     &       + ( 0.334)*r(138)  ! OLE + O3 =
     &       + ( 0.500)*r(141)  ! IOLE + O3 =
     &       + ( 0.120)*r(145)  ! ISO2 + HO2 =
     &       + ( 0.266)*r(149)  ! ISOP + O3 =
     &       + ( 0.461)*r(152)  ! ISPD + O3 =
     &       +          r(156)  ! HPLD =
     &       + ( 1.125)*r(159)  ! EPX2 + HO2 =
     &       + ( 0.125)*r(160)  ! EPX2 + NO =
     &       + ( 0.100)*r(161)  ! EPX2 + C2O3 =
     &       + ( 0.125)*r(162)  ! EPX2 + RO2 =
     &       + ( 0.570)*r(165)  ! TERP + O3 =
     &       + ( 0.500)*r(191)  ! XOPN + O3 =
     &       + ( 0.500)*r(195)  ! OPEN + O3 =
     &       + ( 0.500)*r(202)  ! OPO3 + HO2 =
     &       +          r(210)  ! HOI =
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
     &       + ( 0.170)*r(135)  ! ETH + O3 =
     &       + ( 0.334)*r(138)  ! OLE + O3 =
     &       + ( 0.500)*r(141)  ! IOLE + O3 =
     &       + ( 0.266)*r(149)  ! ISOP + O3 =
     &       + ( 0.461)*r(152)  ! ISPD + O3 =
     &       + ( 0.570)*r(165)  ! TERP + O3 =
     &       + ( 0.500)*r(191)  ! XOPN + O3 =
     &       + ( 0.500)*r(195)  ! OPEN + O3 =
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
     &       +          r( 43)  ! HONO =
     &       +          r( 47)  ! HNO3 =
     &       +          r( 88)  ! MEPX =
     &       +          r( 90)  ! ROOH =
     &       + ( 0.190)*r(111)  ! GLYD =
     &       +          r(156)  ! HPLD =
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
      PA(nn) =          r( 43)  ! HONO =
c
c --- OH from HPLD
c
      nn = nn + 1
      ptname(nn)  = 'newOH_HPLD'
      cpadesc(nn) = 'OH production rate from HPALD photolysis'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =          r(156)  ! HPLD =
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
     &       +          r( 32)  ! NO3 + OH =
     &       +          r( 40)  ! NO + OH =
     &       +          r( 44)  ! HONO + OH =
     &       +          r( 45)  ! NO2 + OH =
     &       +          r( 46)  ! HNO3 + OH =
     &       +          r( 51)  ! PNA + OH =
     &       +          r( 52)  ! SO2 + OH =
     &       + ( 0.600)*r( 87)  ! MEPX + OH =
     &       + ( 0.600)*r( 89)  ! ROOH + OH =
     &       +          r( 91)  ! NTR1 + OH =
     &       +          r( 93)  ! FACD + OH =
     &       +          r( 94)  ! AACD + OH =
     &       +          r( 95)  ! PACD + OH =
     &       +          r( 96)  ! FORM + OH =
     &       +          r(104)  ! ALD2 + OH =
     &       +          r(107)  ! ALDX + OH =
     &       +          r(110)  ! GLYD + OH =
     &       +          r(113)  ! GLY + OH =
     &       +          r(118)  ! MGLY + OH =
     &       +          r(119)  ! H2 + OH =
     &       +          r(120)  ! CO + OH =
     &       +          r(121)  ! CH4 + OH =
     &       +          r(122)  ! ETHA + OH =
     &       +          r(123)  ! MEOH + OH =
     &       +          r(124)  ! ETOH + OH =
     &       +          r(127)  ! ACET + OH =
     &       +          r(128)  ! PRPA + OH =
     &       +          r(129)  ! PAR + OH =
     &       + ( 0.300)*r(133)  ! ETHY + OH =
     &       +          r(134)  ! ETH + OH =
     &       +          r(137)  ! OLE + OH =
     &       +          r(140)  ! IOLE + OH =
     &       +          r(143)  ! ISOP + OH =
     &       +          r(151)  ! ISPD + OH =
     &       + ( 0.067)*r(155)  ! ISPX + OH =
     &       +          r(158)  ! EPOX + OH =
     &       +          r(163)  ! INTR + OH =
     &       +          r(164)  ! TERP + OH =
     &       + ( 0.882)*r(167)  ! BENZ + OH =
     &       + ( 0.900)*r(172)  ! TOL + OH =
     &       + ( 0.756)*r(177)  ! XYL + OH =
     &       +          r(182)  ! CRES + OH =
     &       +          r(186)  ! CRON + OH =
     &       +          r(190)  ! XOPN + OH =
     &       +          r(194)  ! OPEN + OH =
     &       +          r(197)  ! CAT1 + OH =
     &       +          r(205)  ! OPAN + OH =
     &       +          r(206)  ! PANX + OH =
     &       +          r(208)  ! ECH4 + OH =
     &       +          r(218)  ! OIO + OH =
     &       +          r(231)  ! DMS + OH =
     &       +          r(232)  ! DMS + OH + O2 =
     &       +          r(234)  ! OH + NO2 + H2O =
      OH_loss = PA(nn)
c
c --- OH with CO
c
      nn = nn + 1
      ptname(nn)  = 'OHwCO'
      cpadesc(nn) = 'OH reaction rate with CO'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(120)  ! CO + OH =
c
c --- OH with ECH4
c
      nn = nn + 1
      ptname(nn)  = 'OHwECH4'
      cpadesc(nn) = 'OH reaction rate with ECH4 (emitted CH4)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(208)  ! ECH4 + OH =
c
c --- OH with isoprene
c
      nn = nn + 1
      ptname(nn)  = 'OHwISOP'
      cpadesc(nn) = 'OH reaction rate with isoprene'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(143)  ! ISOP + OH =
c
c --- OH with VOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwVOC'
      cpadesc(nn) = 'OH reaction rate with all VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 96)  ! FORM + OH =
     &       +          r(104)  ! ALD2 + OH =
     &       +          r(107)  ! ALDX + OH =
     &       +          r(122)  ! ETHA + OH =
     &       +          r(123)  ! MEOH + OH =
     &       +          r(124)  ! ETOH + OH =
     &       +          r(127)  ! ACET + OH =
     &       +          r(128)  ! PRPA + OH =
     &       +          r(129)  ! PAR + OH =
     &       + ( 0.300)*r(133)  ! ETHY + OH =
     &       +          r(134)  ! ETH + OH =
     &       +          r(137)  ! OLE + OH =
     &       +          r(140)  ! IOLE + OH =
     &       +          r(143)  ! ISOP + OH =
     &       +          r(164)  ! TERP + OH =
     &       + ( 0.882)*r(167)  ! BENZ + OH =
     &       + ( 0.900)*r(172)  ! TOL + OH =
     &       + ( 0.756)*r(177)  ! XYL + OH =
c
c --- OH with HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwHRVOC'
      cpadesc(nn) = 'OH reaction rate with HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(134)  ! ETH + OH =
     &       +          r(137)  ! OLE + OH =
     &       +          r(140)  ! IOLE + OH =
c
c --- OH with Aromatics
c
      nn = nn + 1
      ptname(nn)  = 'OHwArom'
      cpadesc(nn) = 'OH reaction rate with aromatic VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.882)*r(167)  ! BENZ + OH =
     &       + ( 0.900)*r(172)  ! TOL + OH =
     &       + ( 0.756)*r(177)  ! XYL + OH =
c
c --- OH with Alkanes (except methane)
c
      nn = nn + 1
      ptname(nn)  = 'OHwAlkane'
      cpadesc(nn) = 'OH reaction rate with alkanes (except methane)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(122)  ! ETHA + OH =
     &       +          r(128)  ! PRPA + OH =
     &       +          r(129)  ! PAR + OH =
c
c --- Total HCHO production
c
      nn = nn + 1
      ptname(nn)  = 'HCHO_prod'
      cpadesc(nn) = 'Total HCHO production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 71)  ! MEO2 + NO =
     &       + ( 0.100)*r( 72)  ! MEO2 + HO2 =
     &       +          r( 73)  ! MEO2 + C2O3 =
     &       + ( 0.685)*r( 74)  ! MEO2 + RO2 =
     &       + ( 0.400)*r( 87)  ! MEPX + OH =
     &       +          r(101)  ! HCO3 =
     &       + ( 0.740)*r(111)  ! GLYD =
     &       +          r(123)  ! MEOH + OH =
     &       + ( 0.078)*r(124)  ! ETOH + OH =
     &       +          r(127)  ! ACET + OH =
     &       + ( 1.560)*r(134)  ! ETH + OH =
     &       +          r(135)  ! ETH + O3 =
     &       + ( 1.125)*r(136)  ! ETH + NO3 =
     &       + ( 0.781)*r(137)  ! OLE + OH =
     &       + ( 0.555)*r(138)  ! OLE + O3 =
     &       + ( 0.500)*r(139)  ! OLE + NO3 =
     &       + ( 0.128)*r(141)  ! IOLE + O3 =
     &       + ( 0.673)*r(144)  ! ISO2 + NO =
     &       + ( 0.120)*r(145)  ! ISO2 + HO2 =
     &       + ( 0.598)*r(146)  ! ISO2 + C2O3 =
     &       + ( 0.598)*r(147)  ! ISO2 + RO2 =
     &       + ( 0.600)*r(149)  ! ISOP + O3 =
     &       + ( 0.350)*r(150)  ! ISOP + NO3 =
     &       + ( 0.231)*r(152)  ! ISPD + O3 =
     &       + ( 0.260)*r(154)  ! ISPD =
     &       + ( 0.375)*r(159)  ! EPX2 + HO2 =
     &       + ( 0.375)*r(160)  ! EPX2 + NO =
     &       + ( 0.300)*r(161)  ! EPX2 + C2O3 =
     &       + ( 0.375)*r(162)  ! EPX2 + RO2 =
     &       + ( 0.592)*r(163)  ! INTR + OH =
     &       + ( 0.280)*r(164)  ! TERP + OH =
     &       + ( 0.240)*r(165)  ! TERP + O3 =
     &       +          r(188)  ! CRON =
     &       + ( 0.080)*r(195)  ! OPEN + O3 =
     &       + ( 0.140)*r(197)  ! CAT1 + OH =
     &       +          r(231)  ! DMS + OH =
     &       +          r(233)  ! DMS + NO3 =
c
c --- HCHO production from HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'nwHCHO_HRV'
      cpadesc(nn) = 'HCHO production rate from HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 1.560)*r(134)  ! ETH + OH =
     &       +          r(135)  ! ETH + O3 =
     &       + ( 1.125)*r(136)  ! ETH + NO3 =
     &       + ( 0.781)*r(137)  ! OLE + OH =
     &       + ( 0.555)*r(138)  ! OLE + O3 =
     &       + ( 0.500)*r(139)  ! OLE + NO3 =
     &       + ( 0.128)*r(141)  ! IOLE + O3 =
c
c --- HCHO production from Isoprene at first generation
c
      nn = nn + 1
      ptname(nn)  = 'nwHCHO_ISP'
      cpadesc(nn) = 'HCHO production rate from isoprene at first generation'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.600)*r(149)  ! ISOP + O3 =
     &       + ( 0.350)*r(150)  ! ISOP + NO3 =
     &       + ( 0.673)*r(144)  ! ISO2 + NO =
     &       + ( 0.120)*r(145)  ! ISO2 + HO2 =
     &       + ( 0.598)*r(146)  ! ISO2 + C2O3 =
     &       + ( 0.598)*r(147)  ! ISO2 + RO2 =
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
     &       +          r( 32)  ! NO3 + OH =
     &       +          r( 52)  ! SO2 + OH =
     &       +          r( 71)  ! MEO2 + NO =
     &       + ( 0.900)*r( 73)  ! MEO2 + C2O3 =
     &       + ( 0.370)*r( 74)  ! MEO2 + RO2 =
     &       +          r( 75)  ! XO2H + NO =
     &       + ( 0.800)*r( 77)  ! XO2H + C2O3 =
     &       + ( 0.600)*r( 78)  ! XO2H + RO2 =
     &       + ( 0.800)*r( 85)  ! XO2N + C2O3 =
     &       +          r( 90)  ! ROOH =
     &       +          r( 93)  ! FACD + OH =
     &       +          r( 96)  ! FORM + OH =
     &       + ( 2.000)*r( 97)  ! FORM =
     &       +          r( 99)  ! FORM + NO3 =
     &       +          r(101)  ! HCO3 =
     &       +          r(102)  ! HCO3 + NO =
     &       +          r(106)  ! ALD2 =
     &       +          r(109)  ! ALDX =
     &       + ( 0.200)*r(110)  ! GLYD + OH =
     &       + ( 1.400)*r(111)  ! GLYD =
     &       +          r(113)  ! GLY + OH =
     &       + ( 2.000)*r(114)  ! GLY =
     &       +          r(115)  ! GLY + NO3 =
     &       +          r(116)  ! MGLY =
     &       +          r(119)  ! H2 + OH =
     &       +          r(120)  ! CO + OH =
     &       +          r(123)  ! MEOH + OH =
     &       + ( 0.900)*r(124)  ! ETOH + OH =
     &       +          r(131)  ! ROR + O2 =
     &       + ( 0.300)*r(133)  ! ETHY + OH =
     &       + ( 0.270)*r(135)  ! ETH + O3 =
     &       + ( 0.080)*r(138)  ! OLE + O3 =
     &       + ( 0.818)*r(144)  ! ISO2 + NO =
     &       + ( 0.728)*r(146)  ! ISO2 + C2O3 =
     &       + ( 0.728)*r(147)  ! ISO2 + RO2 =
     &       +          r(148)  ! ISO2 =
     &       + ( 0.066)*r(149)  ! ISOP + O3 =
     &       + ( 0.137)*r(151)  ! ISPD + OH =
     &       + ( 0.398)*r(152)  ! ISPD + O3 =
     &       + ( 0.760)*r(154)  ! ISPD =
     &       + ( 0.825)*r(160)  ! EPX2 + NO =
     &       + ( 0.660)*r(161)  ! EPX2 + C2O3 =
     &       + ( 0.825)*r(162)  ! EPX2 + RO2 =
     &       + ( 0.530)*r(167)  ! BENZ + OH =
     &       + ( 0.918)*r(168)  ! BZO2 + NO =
     &       +          r(169)  ! BZO2 + C2O3 =
     &       +          r(171)  ! BZO2 + RO2 =
     &       + ( 0.180)*r(172)  ! TOL + OH =
     &       + ( 0.860)*r(173)  ! TO2 + NO =
     &       +          r(174)  ! TO2 + C2O3 =
     &       +          r(176)  ! TO2 + RO2 =
     &       + ( 0.155)*r(177)  ! XYL + OH =
     &       + ( 0.860)*r(178)  ! XLO2 + NO =
     &       +          r(180)  ! XLO2 + C2O3 =
     &       +          r(181)  ! XLO2 + RO2 =
     &       +          r(182)  ! CRES + OH =
     &       +          r(188)  ! CRON =
     &       + ( 0.700)*r(189)  ! XOPN =
     &       +          r(193)  ! OPEN =
     &       + ( 0.560)*r(195)  ! OPEN + O3 =
     &       + ( 0.200)*r(197)  ! CAT1 + OH =
     &       + ( 0.800)*r(199)  ! OPO3 + NO =
c
c --- HO2 from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_O3'
      cpadesc(nn) = 'HO2 production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.270)*r(135)  ! ETH + O3 =
     &       + ( 0.080)*r(138)  ! OLE + O3 =
     &       + ( 0.066)*r(149)  ! ISOP + O3 =
     &       + ( 0.398)*r(152)  ! ISPD + O3 =
     &       + ( 0.560)*r(195)  ! OPEN + O3 =
      newHO2_O3 = PA(nn)
c
c --- HO2 directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_pht'
      cpadesc(nn) = 'HO2 production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 90)  ! ROOH =
     &       + ( 2.000)*r( 97)  ! FORM =
     &       +          r(101)  ! HCO3 =
     &       +          r(106)  ! ALD2 =
     &       +          r(109)  ! ALDX =
     &       + ( 1.400)*r(111)  ! GLYD =
     &       + ( 2.000)*r(114)  ! GLY =
     &       +          r(116)  ! MGLY =
     &       +          r(148)  ! ISO2 =
     &       + ( 0.760)*r(154)  ! ISPD =
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
      nn = nn + 1
      ptname(nn)  = 'nwHO2_HCHO'
      cpadesc(nn) = 'HO2 production rate from HCHO photolysis'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 97)  ! FORM =
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
     &       +          r(  6)  ! O + NO2 =
     &       +          r( 26)  ! NO2 + O3 =
     &       +          r( 37)  ! N2O5 =
     &       +          r( 38)  ! N2O5 =
     &       +          r( 46)  ! HNO3 + OH =
     &       + ( 0.410)*r( 50)  ! PNA =
     &       + ( 0.400)*r( 56)  ! PAN =
     &       + ( 0.400)*r( 64)  ! PANX =
     &       + ( 0.185)*r(163)  ! INTR + OH =
     &       +          r(223)  ! INO3 =
c
c --- NO3 production from N2O5
c
      nn = nn + 1
      ptname(nn)  = 'N2O5toNO3'
      cpadesc(nn) = 'NO3 production from N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 37)  ! N2O5 =
     &       +          r( 38)  ! N2O5 =
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
     &       +          r( 31)  ! NO3 + O =
     &       +          r( 32)  ! NO3 + OH =
     &       +          r( 33)  ! NO3 + HO2 =
     &       +          r( 34)  ! NO3 + O3 =
     &       + ( 2.000)*r( 35)  ! NO3 + NO3 =
     &       +          r( 36)  ! NO3 + NO2 =
     &       +          r( 99)  ! FORM + NO3 =
     &       +          r(105)  ! ALD2 + NO3 =
     &       +          r(108)  ! ALDX + NO3 =
     &       +          r(112)  ! GLYD + NO3 =
     &       +          r(115)  ! GLY + NO3 =
     &       +          r(117)  ! MGLY + NO3 =
     &       +          r(136)  ! ETH + NO3 =
     &       +          r(139)  ! OLE + NO3 =
     &       +          r(142)  ! IOLE + NO3 =
     &       +          r(150)  ! ISOP + NO3 =
     &       +          r(153)  ! ISPD + NO3 =
     &       +          r(157)  ! HPLD + NO3 =
     &       +          r(166)  ! TERP + NO3 =
     &       +          r(183)  ! CRES + NO3 =
     &       +          r(187)  ! CRON + NO3 =
     &       +          r(192)  ! XOPN + NO3 =
     &       +          r(196)  ! OPEN + NO3 =
     &       +          r(198)  ! CAT1 + NO3 =
     &       +          r(233)  ! DMS + NO3 =
c
c --- NO3 to N2O5
c
      nn = nn + 1
      ptname(nn)  = 'NO3toN2O5'
      cpadesc(nn) = 'NO3 conversion to N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 36)  ! NO3 + NO2 =
c
c --- RO2 loss
c
      nn = nn + 1
      ptname(nn)  = 'RO2_loss'
      cpadesc(nn) = 'RO2 loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 58)  ! C2O3 + RO2 =
     &       +          r( 66)  ! CXO3 + RO2 =
     &       +          r( 68)  ! RO2 + NO =
     &       +          r( 69)  ! RO2 + HO2 =
     &       + ( 2.000)*r( 70)  ! RO2 + RO2 =
      RO2_loss = PA(nn)
c
c --- RO2 with NO
c
      nn = nn + 1
      ptname(nn)  = 'RO2wNO'
      cpadesc(nn) = 'RO2 reaction with NO'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 68)  ! RO2 + NO =
c
c --- RO2 with HO2
c
      nn = nn + 1
      ptname(nn)  = 'RO2wHO2'
      cpadesc(nn) = 'RO2 reaction with HO2'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 69)  ! RO2 + HO2 =
c
c --- RO2 with RO2
c
      nn = nn + 1
      ptname(nn)  = 'RO2wRO2'
      cpadesc(nn) = 'RO2 + RO2 self-reaction'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 70)  ! RO2 + RO2 =
c
c --- Total production of ON compounds
c
      nn = nn + 1
      ptname(nn)  = 'ON_prod'
      cpadesc(nn) = 'Total production of organic nitrate (ON)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.500)*r( 83)  ! XO2N + NO =
     &       +          r(132)  ! ROR + NO2 =
     &       + ( 0.500)*r(136)  ! ETH + NO3 =
     &       + ( 0.500)*r(139)  ! OLE + NO3 =
     &       + ( 0.500)*r(142)  ! IOLE + NO3 =
     &       + ( 0.500)*r( 83)  ! XO2N + NO =
     &       + ( 0.650)*r(150)  ! ISOP + NO3 =
     &       + ( 0.142)*r(153)  ! ISPD + NO3 =
     &       + ( 0.530)*r(166)  ! TERP + NO3 =
     &       + ( 0.082)*r(168)  ! BZO2 + NO =
     &       + ( 0.140)*r(173)  ! TO2 + NO =
     &       + ( 0.140)*r(178)  ! XLO2 + NO =
     &       + ( 0.500)*r(192)  ! XOPN + NO3 =
     &       + ( 0.500)*r(205)  ! OPAN + OH =
     &       + ( 0.100)*r(144)  ! ISO2 + NO =
     &       +          r(184)  ! CRO + NO2 =
      ON_prod = PA(nn)
c
c --- INTR production
c
      nn = nn + 1
      ptname(nn)  = 'INTR_prod'
      cpadesc(nn) = 'Total production of the species INTR'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.100)*r(144)  ! ISO2 + NO =
c
c --- NTR1 production
c
      nn = nn + 1
      ptname(nn)  = 'NTR1_prod'
      cpadesc(nn) = 'Total production of the species NTR1'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.500)*r( 83)  ! XO2N + NO =
     &       +          r(132)  ! ROR + NO2 =
     &       + ( 0.500)*r(136)  ! ETH + NO3 =
     &       + ( 0.500)*r(139)  ! OLE + NO3 =
     &       + ( 0.500)*r(142)  ! IOLE + NO3 =
c
c --- NTR2 production
c
      nn = nn + 1
      ptname(nn)  = 'NTR2_prod'
      cpadesc(nn) = 'Total production of the species NTR2'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.500)*r( 83)  ! XO2N + NO =
     &       +          r( 91)  ! NTR1 + OH =
     &       + ( 0.650)*r(150)  ! ISOP + NO3 =
     &       + ( 0.142)*r(153)  ! ISPD + NO3 =
     &       + ( 0.266)*r(163)  ! INTR + OH =
     &       + ( 0.530)*r(166)  ! TERP + NO3 =
     &       + ( 0.082)*r(168)  ! BZO2 + NO =
     &       + ( 0.140)*r(173)  ! TO2 + NO =
     &       + ( 0.140)*r(178)  ! XLO2 + NO =
     &       +          r(186)  ! CRON + OH =
     &       +          r(187)  ! CRON + NO3 =
     &       + ( 0.500)*r(192)  ! XOPN + NO3 =
     &       + ( 0.500)*r(205)  ! OPAN + OH =
c
c --- NTR2 production from NTR1 with OH
c
      nn = nn + 1
      ptname(nn)  = 'NTR1wOH'
      cpadesc(nn) = 'NTR2 production from NTR1 + OH'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 91)  ! NTR1 + OH =
c
c --- HNO3 production
c
      nn = nn + 1
      ptname(nn)  = 'HNO3_prod'
      cpadesc(nn) = 'HNO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 39)  ! N2O5 + H2O =
     &       +          r( 45)  ! NO2 + OH =
     &       +          r( 99)  ! FORM + NO3 =
     &       +          r(105)  ! ALD2 + NO3 =
     &       +          r(108)  ! ALDX + NO3 =
     &       +          r(112)  ! GLYD + NO3 =
     &       +          r(115)  ! GLY + NO3 =
     &       +          r(117)  ! MGLY + NO3 =
     &       + ( 0.717)*r(153)  ! ISPD + NO3 =
     &       +          r(157)  ! HPLD + NO3 =
     &       +          r(183)  ! CRES + NO3 =
     &       +          r(187)  ! CRON + NO3 =
     &       +          r(196)  ! OPEN + NO3 =
     &       +          r(198)  ! CAT1 + NO3 =
     &       +          r(207)  ! NTR2 =
     &       +          r(224)  ! INO3 + H2O =
     &       +          r(229)  ! INTR =
     &       +          r(233)  ! DMS + NO3 =
     &       +          r(234)  ! OH + NO2 + H2O =
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
     &       + ( 2.000)*r( 39)  ! N2O5 + H2O =
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
     &       +          r( 47)  ! HNO3 =
     &       +          r( 92)  ! NTR1 =
c
c --- Net production of PAN compounds
c
      nn = nn + 1
      ptname(nn)  = 'PAN_prdNet'
      cpadesc(nn) = 'Net production of PAN species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 54)  ! C2O3 + NO2 =
     &       +          r( 62)  ! CXO3 + NO2 =
     &       +          r(200)  ! OPO3 + NO2 =
     &       -          r( 55)  ! PAN =
     &       -          r( 56)  ! PAN =
     &       -          r( 63)  ! PANX =
     &       -          r( 64)  ! PANX =
     &       -          r(201)  ! OPAN =
     &       -          r(205)  ! OPAN + OH =
     &       -          r(206)  ! PANX + OH =
      PAN_prdNet = PA(nn)
c
c --- NO3 with VOC
c
      nn = nn + 1
      ptname(nn)  = 'NO3wVOC'
      cpadesc(nn) = 'NO3 reaction with VOC species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 99)  ! FORM + NO3 =
     &       +          r(105)  ! ALD2 + NO3 =
     &       +          r(108)  ! ALDX + NO3 =
     &       +          r(136)  ! ETH + NO3 =
     &       +          r(139)  ! OLE + NO3 =
     &       +          r(142)  ! IOLE + NO3 =
     &       +          r(150)  ! ISOP + NO3 =
     &       +          r(166)  ! TERP + NO3 =
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
     &       + ( 2.000)*r(209)  ! I2 =
     &       +          r(210)  ! HOI =
     &       +          r(212)  ! IO =
     &       + ( 0.400)*r(213)  ! IO + IO =
     &       +          r(215)  ! IO + NO =
     &       +          r(217)  ! OIO =
     &       +          r(221)  ! I2O2 =
     &       +          r(223)  ! INO3 =
c
c --- Ozone destruction by I
c
      nn = nn + 1
      ptname(nn)  = 'I_O3dest'
      cpadesc(nn) = 'Ozone destruction by I atom'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       -          r(214)  ! IO + HO2 =
     &       - ( 0.400)*r(213)  ! IO + IO =
     &       -          r(216)  ! IO + NO2 =
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

