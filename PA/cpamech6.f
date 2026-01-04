      subroutine cpamech6(r, rk, dtfact, nr, pa, npa, npa_init, ldark )

      use filunit
      use tracer
      use procan
      implicit none
c
c----CAMx v7.20 220430
c
c     CPAMECH6 computes chemical process analysis (CPA) parameters
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
      integer  n, nn, light
      real     sum, po3, do3, HO2toNO2, HO2_loss
      real     NO2wOH, HO2wHO2, ratio_ind
      real     OH_new, HO2_new, OH_loss, RO2_loss, HOx_CL
      real     newOH_O1D, newOH_O3, newOH_phot
      real     newHO2_O3, newHO2_pht
      real     NO3wVOC, N2O5wH2O, PAN_prdNet, ON_prod
c
c --- Entry point:
c
      nn = 0
      light = 1
      if (ldark) light = 0
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
     &       + ( 0.200)*r( 92)  ! C2O3 + HO2 =
     &       + ( 0.200)*r(108)  ! CXO3 + HO2 =
     &       -          r(  3)  ! O3 + NO =
     &       -          r(  7)  ! NO2 + O3 =
     &       -          r(  8)  ! O3 =
     &       -          r(  9)  ! O3 =
     &       -          r( 12)  ! O3 + OH =
     &       -          r( 13)  ! O3 + HO2 =
     &       -          r( 49)  ! NO3 + O3 =
     &       -          r(121)  ! O3 + OLE =
     &       -          r(125)  ! O3 + ETH =
     &       -          r(129)  ! IOLE + O3 =
     &       -          r(140)  ! OPEN + O3 =
     &       -          r(146)  ! O3 + ISOP =
     &       -          r(150)  ! O3 + ISPD =
     &       -          r(155)  ! TERP + O3 =
c
      po3 = max(0.0, PA(nn))
      do3 = min(0.0, PA(nn))
c
c --- Calculate the P(H2O2)/P(HNO3) indicator ratio and apply to PO3
c
      HO2wHO2 = + r( 34) + r( 35)
      NO2wOH = + r( 28)
      ratio_ind = min( 10., HO2wHO2/max( 1.0E-12, NO2wOH ) )*light
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
     &       +          r( 30)  ! HO2 + NO =
c
      HO2_loss = 
     &       +          r( 13)  ! O3 + HO2 =
     &       +          r( 30)  ! HO2 + NO =
     &       + ( 2.000)*r( 34)  ! HO2 + HO2 =
     &       + ( 2.000)*r( 35)  ! HO2 + HO2 + H2O =
     &       +          r( 43)  ! OH + HO2 =
     &       +          r( 44)  ! HO2 + O =
     &       +          r( 48)  ! NO3 + HO2 =
     &       +          r( 56)  ! XO2 + HO2 =
     &       +          r( 57)  ! XO2N + HO2 =
     &       +          r( 69)  ! MEO2 + HO2 =
     &       +          r( 79)  ! FORM + HO2 =
     &       +          r( 82)  ! HCO3 + HO2 =
     &       +          r( 92)  ! C2O3 + HO2 =
     &       +          r(108)  ! CXO3 + HO2 =
     &       +          r(137)  ! CRO + HO2 =
c
      PA(nn) = PA(nn) - (
     &       + r( 12)  ! O3 + OH =
     &                     ) * (HO2_loss-HO2toNO2)/(HO2_loss+10E-12)
c
c  +  O3 + VOC
c
      PA(nn) = PA(nn)
     &       - r(121)  ! O3 + OLE =
     &       - r(125)  ! O3 + ETH =
     &       - r(129)  ! IOLE + O3 =
     &       - r(140)  ! OPEN + O3 =
     &       - r(146)  ! O3 + ISOP =
     &       - r(150)  ! O3 + ISPD =
     &       - r(155)  ! TERP + O3 =
c
c  +  O(3P) + VOC
c
      PA(nn) = PA(nn)
     &       - r( 77)  ! FORM + O =
     &       - r( 84)  ! ALD2 + O =
     &       - r( 99)  ! ALDX + O =
     &       - r(119)  ! O + OLE =
     &       - r(123)  ! O + ETH =
     &       - r(127)  ! IOLE + O =
     &       - r(144)  ! O + ISOP =
     &       - r(153)  ! TERP + O =
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
     &       +          r( 25)  ! HONO =
     &       +          r( 30)  ! HO2 + NO =
     &       + ( 2.000)*r( 36)  ! H2O2 =
     &       +          r( 38)  ! O1D + H2 =
     &       +          r( 44)  ! HO2 + O =
     &       +          r( 45)  ! H2O2 + O =
     &       + ( 0.390)*r( 51)  ! PNA =
     &       +          r( 52)  ! HNO3 =
     &       +          r( 65)  ! ROOH =
     &       +          r( 72)  ! MEPX =
     &       +          r( 77)  ! FORM + O =
     &       +          r( 84)  ! ALD2 + O =
     &       +          r( 97)  ! PACD =
     &       +          r( 99)  ! ALDX + O =
     &       + ( 0.100)*r(119)  ! O + OLE =
     &       + ( 0.100)*r(121)  ! O3 + OLE =
     &       + ( 0.300)*r(123)  ! O + ETH =
     &       + ( 0.130)*r(125)  ! O3 + ETH =
     &       + ( 0.500)*r(129)  ! IOLE + O3 =
     &       + ( 0.080)*r(140)  ! OPEN + O3 =
     &       + ( 0.266)*r(146)  ! O3 + ISOP =
     &       + ( 0.268)*r(150)  ! O3 + ISPD =
     &       + ( 0.570)*r(155)  ! TERP + O3 =
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
     &       + ( 0.100)*r(121)  ! O3 + OLE =
     &       + ( 0.130)*r(125)  ! O3 + ETH =
     &       + ( 0.500)*r(129)  ! IOLE + O3 =
     &       + ( 0.080)*r(140)  ! OPEN + O3 =
     &       + ( 0.266)*r(146)  ! O3 + ISOP =
     &       + ( 0.268)*r(150)  ! O3 + ISPD =
     &       + ( 0.570)*r(155)  ! TERP + O3 =
      newOH_O3 = PA(nn)
c
c --- OH directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newOH_phot'
      cpadesc(nn) = 'OH production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 25)  ! HONO =
     &       + ( 2.000)*r( 36)  ! H2O2 =
     &       +          r( 52)  ! HNO3 =
     &       +          r( 65)  ! ROOH =
     &       +          r( 72)  ! MEPX =
     &       +          r( 97)  ! PACD =
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
      PA(nn) =          r( 25)  ! HONO =
c
c --- OH loss
c
      nn = nn + 1
      ptname(nn)  = 'OH_loss'
      cpadesc(nn) = 'OH loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 12)  ! O3 + OH =
     &       +          r( 24)  ! NO + OH =
     &       +          r( 26)  ! OH + HONO =
     &       +          r( 28)  ! NO2 + OH =
     &       +          r( 29)  ! OH + HNO3 =
     &       +          r( 33)  ! OH + PNA =
     &       +          r( 37)  ! OH + H2O2 =
     &       +          r( 39)  ! OH + H2 =
     &       +          r( 40)  ! OH + O =
     &       + ( 2.000)*r( 41)  ! OH + OH =
     &       + ( 2.000)*r( 42)  ! OH + OH =
     &       +          r( 43)  ! OH + HO2 =
     &       +          r( 47)  ! NO3 + OH =
     &       +          r( 61)  ! NTR + OH =
     &       +          r( 63)  ! SO2 + OH =
     &       +          r( 64)  ! ROOH + OH =
     &       +          r( 66)  ! OH + CO =
     &       +          r( 67)  ! OH + CH4 =
     &       +          r( 71)  ! MEPX + OH =
     &       +          r( 73)  ! MEOH + OH =
     &       +          r( 74)  ! FORM + OH =
     &       +          r( 83)  ! FACD + OH =
     &       +          r( 85)  ! ALD2 + OH =
     &       +          r( 96)  ! PACD + OH =
     &       +          r( 98)  ! AACD + OH =
     &       +          r(100)  ! ALDX + OH =
     &       +          r(107)  ! PANX + OH =
     &       +          r(113)  ! OH + ETHA =
     &       +          r(114)  ! OH + ETOH =
     &       +          r(115)  ! PAR + OH =
     &       +          r(120)  ! OH + OLE =
     &       +          r(124)  ! OH + ETH =
     &       +          r(128)  ! IOLE + OH =
     &       +          r(131)  ! TOL + OH =
     &       +          r(134)  ! OH + CRES =
     &       +          r(139)  ! OPEN + OH =
     &       +          r(141)  ! OH + XYL =
     &       +          r(142)  ! OH + MGLY =
     &       +          r(145)  ! OH + ISOP =
     &       +          r(149)  ! OH + ISPD =
     &       +          r(154)  ! TERP + OH =
      OH_loss = PA(nn)
c
c --- OH with CO
c
      nn = nn + 1
      ptname(nn)  = 'OHwCO'
      cpadesc(nn) = 'OH reaction rate with CO'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 66)  ! OH + CO =
c
c --- OH with isoprene
c
      nn = nn + 1
      ptname(nn)  = 'OHwISOP'
      cpadesc(nn) = 'OH reaction rate with isoprene'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(145)  ! OH + ISOP =
c
c --- OH with VOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwVOC'
      cpadesc(nn) = 'OH reaction rate with all VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 73)  ! MEOH + OH =
     &       +          r( 74)  ! FORM + OH =
     &       +          r( 85)  ! ALD2 + OH =
     &       +          r(100)  ! ALDX + OH =
     &       +          r(113)  ! OH + ETHA =
     &       +          r(114)  ! OH + ETOH =
     &       +          r(115)  ! PAR + OH =
     &       +          r(120)  ! OH + OLE =
     &       +          r(124)  ! OH + ETH =
     &       +          r(128)  ! IOLE + OH =
     &       +          r(131)  ! TOL + OH =
     &       +          r(141)  ! OH + XYL =
     &       +          r(145)  ! OH + ISOP =
     &       +          r(154)  ! TERP + OH =
c
c --- OH with HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwHRVOC'
      cpadesc(nn) = 'OH reaction rate with HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(120)  ! OH + OLE =
     &       +          r(124)  ! OH + ETH =
     &       +          r(128)  ! IOLE + OH =
c
c --- OH with Aromatics
c
      nn = nn + 1
      ptname(nn)  = 'OHwArom'
      cpadesc(nn) = 'OH reaction rate with aromatic VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(131)  ! TOL + OH =
     &       +          r(141)  ! OH + XYL =
c
c --- OH with Alkanes (except methane)
c
      nn = nn + 1
      ptname(nn)  = 'OHwAlkane'
      cpadesc(nn) = 'OH reaction rate with alkanes (except methane)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(113)  ! OH + ETHA =
     &       +          r(115)  ! PAR + OH =
c
c --- Total HCHO production
c
      nn = nn + 1
      ptname(nn)  = 'HCHO_prod'
      cpadesc(nn) = 'Total HCHO production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.330)*r( 61)  ! NTR + OH =
     &       + ( 0.330)*r( 62)  ! NTR =
     &       +          r( 68)  ! MEO2 + NO =
     &       + ( 1.370)*r( 70)  ! MEO2 + MEO2 =
     &       +          r( 72)  ! MEPX =
     &       +          r( 73)  ! MEOH + OH =
     &       +          r( 80)  ! HCO3 =
     &       +          r( 93)  ! C2O3 + MEO2 =
     &       + ( 0.100)*r(109)  ! CXO3 + MEO2 =
     &       + ( 0.100)*r(114)  ! OH + ETOH =
     &       + ( 0.200)*r(119)  ! O + OLE =
     &       + ( 0.800)*r(120)  ! OH + OLE =
     &       + ( 0.740)*r(121)  ! O3 + OLE =
     &       +          r(122)  ! NO3 + OLE =
     &       +          r(123)  ! O + ETH =
     &       + ( 1.560)*r(124)  ! OH + ETH =
     &       +          r(125)  ! O3 + ETH =
     &       + ( 2.000)*r(126)  ! NO3 + ETH =
     &       + ( 0.250)*r(129)  ! IOLE + O3 =
     &       +          r(139)  ! OPEN + OH =
     &       + ( 0.700)*r(140)  ! OPEN + O3 =
     &       + ( 0.500)*r(144)  ! O + ISOP =
     &       + ( 0.629)*r(145)  ! OH + ISOP =
     &       + ( 0.600)*r(146)  ! O3 + ISOP =
     &       + ( 0.167)*r(149)  ! OH + ISPD =
     &       + ( 0.150)*r(150)  ! O3 + ISPD =
     &       + ( 0.282)*r(151)  ! NO3 + ISPD =
     &       + ( 0.900)*r(152)  ! ISPD =
     &       + ( 0.280)*r(154)  ! TERP + OH =
     &       + ( 0.240)*r(155)  ! TERP + O3 =
c
c --- HCHO production from HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'nwHCHO_HRV'
      cpadesc(nn) = 'HCHO production rate from HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.200)*r(119)  ! O + OLE =
     &       + ( 0.800)*r(120)  ! OH + OLE =
     &       + ( 0.740)*r(121)  ! O3 + OLE =
     &       +          r(122)  ! NO3 + OLE =
     &       +          r(123)  ! O + ETH =
     &       + ( 1.560)*r(124)  ! OH + ETH =
     &       +          r(125)  ! O3 + ETH =
     &       + ( 2.000)*r(126)  ! NO3 + ETH =
     &       + ( 0.250)*r(129)  ! IOLE + O3 =
c
c --- HCHO production from Isoprene at first generation
c
      nn = nn + 1
      ptname(nn)  = 'nwHCHO_ISP'
      cpadesc(nn) = 'HCHO production rate from isoprene at first generation'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.500)*r(144)  ! O + ISOP =
     &       + ( 0.629)*r(145)  ! OH + ISOP =
     &       + ( 0.600)*r(146)  ! O3 + ISOP =
c
c --- Total HO2 production, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'HO2_prod'
      cpadesc(nn) = 'Total HO2 production rate, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 12)  ! O3 + OH =
     &       +          r( 37)  ! OH + H2O2 =
     &       +          r( 38)  ! O1D + H2 =
     &       +          r( 39)  ! OH + H2 =
     &       +          r( 40)  ! OH + O =
     &       +          r( 45)  ! H2O2 + O =
     &       +          r( 47)  ! NO3 + OH =
     &       +          r( 61)  ! NTR + OH =
     &       +          r( 62)  ! NTR =
     &       +          r( 63)  ! SO2 + OH =
     &       +          r( 65)  ! ROOH =
     &       +          r( 66)  ! OH + CO =
     &       +          r( 68)  ! MEO2 + NO =
     &       + ( 0.740)*r( 70)  ! MEO2 + MEO2 =
     &       + ( 0.300)*r( 71)  ! MEPX + OH =
     &       +          r( 72)  ! MEPX =
     &       +          r( 73)  ! MEOH + OH =
     &       +          r( 74)  ! FORM + OH =
     &       + ( 2.000)*r( 75)  ! FORM =
     &       +          r( 77)  ! FORM + O =
     &       +          r( 78)  ! FORM + NO3 =
     &       +          r( 80)  ! HCO3 =
     &       +          r( 81)  ! HCO3 + NO =
     &       +          r( 83)  ! FACD + OH =
     &       +          r( 87)  ! ALD2 =
     &       + ( 0.900)*r( 93)  ! C2O3 + MEO2 =
     &       +          r(102)  ! ALDX =
     &       +          r(103)  ! CXO3 + NO =
     &       +          r(109)  ! CXO3 + MEO2 =
     &       + ( 2.000)*r(111)  ! CXO3 + CXO3 =
     &       +          r(112)  ! CXO3 + C2O3 =
     &       +          r(113)  ! OH + ETHA =
     &       +          r(114)  ! OH + ETOH =
     &       + ( 0.110)*r(115)  ! PAR + OH =
     &       + ( 0.940)*r(116)  ! ROR =
     &       +          r(117)  ! ROR =
     &       + ( 0.300)*r(119)  ! O + OLE =
     &       + ( 0.950)*r(120)  ! OH + OLE =
     &       + ( 0.440)*r(121)  ! O3 + OLE =
     &       + ( 1.700)*r(123)  ! O + ETH =
     &       +          r(124)  ! OH + ETH =
     &       + ( 0.130)*r(125)  ! O3 + ETH =
     &       + ( 0.100)*r(127)  ! IOLE + O =
     &       +          r(128)  ! IOLE + OH =
     &       + ( 0.500)*r(129)  ! IOLE + O3 =
     &       +          r(130)  ! IOLE + NO3 =
     &       + ( 0.440)*r(131)  ! TOL + OH =
     &       + ( 0.900)*r(132)  ! TO2 + NO =
     &       +          r(133)  ! TO2 =
     &       + ( 0.600)*r(134)  ! OH + CRES =
     &       +          r(138)  ! OPEN =
     &       + ( 2.000)*r(139)  ! OPEN + OH =
     &       + ( 0.760)*r(140)  ! OPEN + O3 =
     &       + ( 0.700)*r(141)  ! OH + XYL =
     &       +          r(143)  ! MGLY =
     &       + ( 0.250)*r(144)  ! O + ISOP =
     &       + ( 0.912)*r(145)  ! OH + ISOP =
     &       + ( 0.066)*r(146)  ! O3 + ISOP =
     &       + ( 0.800)*r(147)  ! NO3 + ISOP =
     &       + ( 0.800)*r(148)  ! NO2 + ISOP =
     &       + ( 0.503)*r(149)  ! OH + ISPD =
     &       + ( 0.154)*r(150)  ! O3 + ISPD =
     &       + ( 0.925)*r(151)  ! NO3 + ISPD =
     &       + ( 1.033)*r(152)  ! ISPD =
     &       + ( 0.750)*r(154)  ! TERP + OH =
     &       + ( 0.070)*r(155)  ! TERP + O3 =
     &       + ( 0.280)*r(156)  ! TERP + NO3 =
c
c --- HO2 from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_O3'
      cpadesc(nn) = 'HO2 production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.440)*r(121)  ! O3 + OLE =
     &       + ( 0.130)*r(125)  ! O3 + ETH =
     &       + ( 0.500)*r(129)  ! IOLE + O3 =
     &       + ( 0.760)*r(140)  ! OPEN + O3 =
     &       + ( 0.066)*r(146)  ! O3 + ISOP =
     &       + ( 0.154)*r(150)  ! O3 + ISPD =
     &       + ( 0.070)*r(155)  ! TERP + O3 =
      newHO2_O3 = PA(nn)
c
c --- HO2 directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_pht'
      cpadesc(nn) = 'HO2 production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 62)  ! NTR =
     &       +          r( 65)  ! ROOH =
     &       +          r( 72)  ! MEPX =
     &       + ( 2.000)*r( 75)  ! FORM =
     &       +          r( 87)  ! ALD2 =
     &       +          r(102)  ! ALDX =
     &       + ( 0.940)*r(116)  ! ROR =
     &       +          r(117)  ! ROR =
     &       +          r(133)  ! TO2 =
     &       +          r(138)  ! OPEN =
     &       +          r(143)  ! MGLY =
     &       + ( 1.033)*r(152)  ! ISPD =
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
     &       + ( 2.000)*r( 75)  ! FORM =
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
     &           *light*(dtfact/ppbfact)
      endif
c
c --- NO3 production
c
      nn = nn + 1
      ptname(nn)  = 'NO3_prod'
      cpadesc(nn) = 'NO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  5)  ! O + NO2 =
     &       +          r(  7)  ! NO2 + O3 =
     &       +          r( 21)  ! N2O5 =
     &       +          r( 29)  ! OH + HNO3 =
     &       + ( 0.390)*r( 51)  ! PNA =
     &       +          r( 53)  ! N2O5 =
c
c --- NO3 production from N2O5
c
      nn = nn + 1
      ptname(nn)  = 'N2O5toNO3'
      cpadesc(nn) = 'NO3 production from N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 21)  ! N2O5 =
     &       +          r( 53)  ! N2O5 =
c
c --- NO3 loss
c
      nn = nn + 1
      ptname(nn)  = 'NO3_loss'
      cpadesc(nn) = 'NO3 loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 14)  ! NO3 =
     &       +          r( 15)  ! NO3 =
     &       +          r( 16)  ! NO3 + NO =
     &       +          r( 17)  ! NO3 + NO2 =
     &       +          r( 18)  ! NO3 + NO2 =
     &       +          r( 46)  ! NO3 + O =
     &       +          r( 47)  ! NO3 + OH =
     &       +          r( 48)  ! NO3 + HO2 =
     &       +          r( 49)  ! NO3 + O3 =
     &       + ( 2.000)*r( 50)  ! NO3 + NO3 =
     &       +          r( 78)  ! FORM + NO3 =
     &       +          r( 86)  ! ALD2 + NO3 =
     &       +          r(101)  ! ALDX + NO3 =
     &       +          r(122)  ! NO3 + OLE =
     &       +          r(126)  ! NO3 + ETH =
     &       +          r(130)  ! IOLE + NO3 =
     &       +          r(135)  ! CRES + NO3 =
     &       +          r(147)  ! NO3 + ISOP =
     &       +          r(151)  ! NO3 + ISPD =
     &       +          r(156)  ! TERP + NO3 =
c
c --- NO3 to N2O5
c
      nn = nn + 1
      ptname(nn)  = 'NO3toN2O5'
      cpadesc(nn) = 'NO3 conversion to N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 18)  ! NO3 + NO2 =
c
c --- Total production of ON compounds
c
      nn = nn + 1
      ptname(nn)  = 'ON_prod'
      cpadesc(nn) = 'Total production of organic nitrate (ON)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 55)  ! XO2N + NO =
     &       +          r(118)  ! ROR + NO2 =
     &       + ( 0.100)*r(132)  ! TO2 + NO =
     &       +          r(136)  ! CRO + NO2 =
     &       + ( 0.800)*r(147)  ! NO3 + ISOP =
     &       + ( 0.800)*r(148)  ! NO2 + ISOP =
     &       + ( 0.850)*r(151)  ! NO3 + ISPD =
     &       + ( 0.530)*r(156)  ! TERP + NO3 =
      ON_prod = PA(nn)
c
c --- NTR production
c
      nn = nn + 1
      ptname(nn)  = 'NTR_prod'
      cpadesc(nn) = 'Total production of the species NTR'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 55)  ! XO2N + NO =
     &       +          r(118)  ! ROR + NO2 =
     &       + ( 0.100)*r(132)  ! TO2 + NO =
     &       +          r(136)  ! CRO + NO2 =
     &       + ( 0.800)*r(147)  ! NO3 + ISOP =
     &       + ( 0.800)*r(148)  ! NO2 + ISOP =
     &       + ( 0.850)*r(151)  ! NO3 + ISPD =
     &       + ( 0.530)*r(156)  ! TERP + NO3 =
c
c --- HNO3 production
c
      nn = nn + 1
      ptname(nn)  = 'HNO3_prod'
      cpadesc(nn) = 'HNO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 19)  ! N2O5 + H2O =
     &       + ( 2.000)*r( 20)  ! N2O5 + H2O + H2O =
     &       +          r( 28)  ! NO2 + OH =
     &       +          r( 48)  ! NO3 + HO2 =
     &       +          r( 61)  ! NTR + OH =
     &       +          r( 78)  ! FORM + NO3 =
     &       +          r( 86)  ! ALD2 + NO3 =
     &       +          r(101)  ! ALDX + NO3 =
     &       +          r(135)  ! CRES + NO3 =
     &       + ( 0.150)*r(151)  ! NO3 + ISPD =
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
     &       + ( 2.000)*r( 19)  ! N2O5 + H2O =
     &       + ( 2.000)*r( 20)  ! N2O5 + H2O + H2O =
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
     &       +          r( 52)  ! HNO3 =
     &       +          r( 62)  ! NTR =
c
c --- Net production of PAN compounds
c
      nn = nn + 1
      ptname(nn)  = 'PAN_prdNet'
      cpadesc(nn) = 'Net production of PAN species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 89)  ! C2O3 + NO2 =
     &       +          r(104)  ! CXO3 + NO2 =
     &       -          r( 90)  ! PAN =
     &       -          r( 91)  ! PAN =
     &       -          r(105)  ! PANX =
     &       -          r(106)  ! PANX =
     &       -          r(107)  ! PANX + OH =
      PAN_prdNet = PA(nn)
c
c --- NO3 with VOC
c
      nn = nn + 1
      ptname(nn)  = 'NO3wVOC'
      cpadesc(nn) = 'NO3 reaction with VOC species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 78)  ! FORM + NO3 =
     &       +          r( 86)  ! ALD2 + NO3 =
     &       +          r(101)  ! ALDX + NO3 =
     &       +          r(122)  ! NO3 + OLE =
     &       +          r(126)  ! NO3 + ETH =
     &       +          r(130)  ! IOLE + NO3 =
     &       +          r(147)  ! NO3 + ISOP =
     &       +          r(156)  ! TERP + NO3 =
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
        PA(nn) = max( 0., min( 10., PA(nn) ))*light*(dtfact/ppbfact)
      endif
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

