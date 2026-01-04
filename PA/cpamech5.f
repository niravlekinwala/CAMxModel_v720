      subroutine cpamech5(r, rk, dtfact, nr, pa, npa, npa_init, ldark )

      use filunit
      use tracer
      use procan
      implicit none
c
c----CAMx v7.20 220430
c
c     CPAMECH5 computes chemical process analysis (CPA) parameters
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
      PA(nn) =      rk( 18)*dtfact
c
c
c --- Net O3 production
c
      nn = nn + 1
      ptname(nn)  = 'PO3_net'
      cpadesc(nn) = 'Net O3 production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  2)  ! O3P + O2 + M =
     &       + ( 0.300)*r( 67)  ! MCO3 + HO2 =
     &       + ( 0.250)*r( 77)  ! RCO3 + HO2 =
     &       + ( 0.250)*r( 88)  ! BZC3 + HO2 =
     &       + ( 0.250)*r(100)  ! MAC3 + HO2 =
     &       -          r(  3)  ! O3P + O3 =
     &       -          r(  7)  ! O3 + NO =
     &       -          r(  8)  ! O3 + NO2 =
     &       -          r( 18)  ! O3 =
     &       -          r( 19)  ! O3 =
     &       -          r( 30)  ! OH + O3 =
     &       -          r( 36)  ! HO2 + O3 =
     &       -          r(247)  ! AFG1 + O3 =
     &       -          r(250)  ! AFG2 + O3 =
     &       -          r(253)  ! AFG3 + O3 =
     &       -          r(255)  ! MACR + O3 =
     &       -          r(260)  ! MVK + O3 =
     &       -          r(264)  ! IPRD + O3 =
     &       -          r(275)  ! ACRO + O3 =
     &       -          r(515)  ! ETHE + O3 =
     &       -          r(519)  ! PRPE + O3 =
     &       -          r(523)  ! BD13 + O3 =
     &       -          r(527)  ! ISOP + O3 =
     &       -          r(531)  ! APIN + O3 =
     &       -          r(535)  ! ACYE + O3 =
     &       -          r(549)  ! OLE1 + O3 =
     &       -          r(553)  ! OLE2 + O3 =
     &       -          r(559)  ! TERP + O3 =
     &       -          r(563)  ! SESQ + O3 =
c
      po3 = max(0.0, PA(nn))
      do3 = min(0.0, PA(nn))
c
c --- Calculate the P(H2O2)/P(HNO3) indicator ratio and apply to PO3
c
      HO2wHO2 = + r( 37) + r( 38)
      NO2wOH = + r( 25)
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
     &       - r( 20)  ! O1D + H2O =
c
c  +  HO2 + O3 (assuming no OH recycled)
c
      PA(nn) = PA(nn)
     &       - r( 36)  ! HO2 + O3 =
c
c  +  OH + O3 (accounting for HO2 recycled to O3 via NO2 produced)
c
      HO2toNO2 = 
     &       +          r( 31)  ! HO2 + NO =
     &       + ( 0.800)*r( 39)  ! NO3 + HO2 =
c
      HO2_loss = 
     &       +          r( 31)  ! HO2 + NO =
     &       +          r( 36)  ! HO2 + O3 =
     &       + ( 2.000)*r( 37)  ! HO2 + HO2 =
     &       + ( 2.000)*r( 38)  ! HO2 + HO2 + H2O =
     &       +          r( 39)  ! NO3 + HO2 =
     &       +          r( 43)  ! OH + HO2 =
     &       +          r( 47)  ! MEO2 + HO2 =
     &       +          r( 48)  ! MEO2 + HO2 =
     &       +          r( 53)  ! RO2C + HO2 =
     &       +          r( 58)  ! RO2X + HO2 =
     &       +          r( 67)  ! MCO3 + HO2 =
     &       +          r( 77)  ! RCO3 + HO2 =
     &       +          r( 88)  ! BZC3 + HO2 =
     &       +          r(100)  ! MAC3 + HO2 =
     &       +          r(112)  ! BZO + HO2 =
c
      PA(nn) = PA(nn) - (
     &       + r( 30)  ! OH + O3 =
     &                     ) * (HO2_loss-HO2toNO2)/(HO2_loss+10E-12)
c
c  +  O3 + VOC
c
      PA(nn) = PA(nn)
     &       - r(247)  ! AFG1 + O3 =
     &       - r(250)  ! AFG2 + O3 =
     &       - r(253)  ! AFG3 + O3 =
     &       - r(255)  ! MACR + O3 =
     &       - r(260)  ! MVK + O3 =
     &       - r(264)  ! IPRD + O3 =
     &       - r(275)  ! ACRO + O3 =
     &       - r(515)  ! ETHE + O3 =
     &       - r(519)  ! PRPE + O3 =
     &       - r(523)  ! BD13 + O3 =
     &       - r(527)  ! ISOP + O3 =
     &       - r(531)  ! APIN + O3 =
     &       - r(535)  ! ACYE + O3 =
     &       - r(549)  ! OLE1 + O3 =
     &       - r(553)  ! OLE2 + O3 =
     &       - r(559)  ! TERP + O3 =
     &       - r(563)  ! SESQ + O3 =
c
c  +  O(3P) + VOC
c
      PA(nn) = PA(nn)
     &       - r(257)  ! MACR + O3P =
     &       - r(261)  ! MVK + O3P =
     &       - r(277)  ! ACRO + O3P =
     &       - r(517)  ! ETHE + O3P =
     &       - r(521)  ! PRPE + O3P =
     &       - r(525)  ! BD13 + O3P =
     &       - r(529)  ! ISOP + O3P =
     &       - r(533)  ! APIN + O3P =
     &       - r(551)  ! OLE1 + O3P =
     &       - r(555)  ! OLE2 + O3P =
     &       - r(561)  ! TERP + O3P =
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
     &       + ( 2.000)*r( 20)  ! O1D + H2O =
     &       +          r( 23)  ! HONO =
     &       +          r( 28)  ! HNO3 =
     &       +          r( 31)  ! HO2 + NO =
     &       + ( 0.390)*r( 34)  ! PNA =
     &       +          r( 36)  ! HO2 + O3 =
     &       + ( 0.800)*r( 39)  ! NO3 + HO2 =
     &       + ( 2.000)*r( 41)  ! H2O2 =
     &       +          r(124)  ! XOH + NO =
     &       +          r(126)  ! XOH + NO3 =
     &       + ( 0.500)*r(127)  ! XOH + MEO2 =
     &       + ( 0.500)*r(128)  ! XOH + RO2C =
     &       + ( 0.500)*r(129)  ! XOH + RO2X =
     &       +          r(130)  ! XOH + MCO3 =
     &       +          r(131)  ! XOH + RCO3 =
     &       +          r(132)  ! XOH + BZC3 =
     &       +          r(133)  ! XOH + MAC3 =
     &       +          r(223)  ! COOH =
     &       +          r(225)  ! ROOH =
     &       +          r(227)  ! R6PX =
     &       +          r(229)  ! RAPX =
     &       + ( 0.826)*r(247)  ! AFG1 + O3 =
     &       + ( 0.826)*r(250)  ! AFG2 + O3 =
     &       + ( 0.471)*r(253)  ! AFG3 + O3 =
     &       + ( 0.208)*r(255)  ! MACR + O3 =
     &       + ( 0.330)*r(258)  ! MACR =
     &       + ( 0.164)*r(260)  ! MVK + O3 =
     &       + ( 0.285)*r(264)  ! IPRD + O3 =
     &       + ( 0.330)*r(275)  ! ACRO + O3 =
     &       + ( 0.178)*r(278)  ! ACRO =
     &       +          r(280)  ! CO3H =
     &       +          r(282)  ! RO3H =
     &       + ( 0.160)*r(515)  ! ETHE + O3 =
     &       + ( 0.350)*r(519)  ! PRPE + O3 =
     &       + ( 0.080)*r(523)  ! BD13 + O3 =
     &       + ( 0.266)*r(527)  ! ISOP + O3 =
     &       + ( 0.728)*r(531)  ! APIN + O3 =
     &       + ( 0.500)*r(535)  ! ACYE + O3 =
     &       + ( 0.128)*r(549)  ! OLE1 + O3 =
     &       + ( 0.443)*r(553)  ! OLE2 + O3 =
     &       + ( 0.499)*r(559)  ! TERP + O3 =
     &       + ( 0.499)*r(563)  ! SESQ + O3 =
c
c --- OH from O(1D)
c
      nn = nn + 1
      ptname(nn)  = 'newOH_O1D'
      cpadesc(nn) = 'OH production rate from O(1D) + H2O'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) = ( 2.000)*r( 20)  ! O1D + H2O =
      newOH_O1D = PA(nn)
c
c --- OH from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newOH_O3'
      cpadesc(nn) = 'OH production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.826)*r(247)  ! AFG1 + O3 =
     &       + ( 0.826)*r(250)  ! AFG2 + O3 =
     &       + ( 0.471)*r(253)  ! AFG3 + O3 =
     &       + ( 0.208)*r(255)  ! MACR + O3 =
     &       + ( 0.164)*r(260)  ! MVK + O3 =
     &       + ( 0.285)*r(264)  ! IPRD + O3 =
     &       + ( 0.330)*r(275)  ! ACRO + O3 =
     &       + ( 0.160)*r(515)  ! ETHE + O3 =
     &       + ( 0.350)*r(519)  ! PRPE + O3 =
     &       + ( 0.080)*r(523)  ! BD13 + O3 =
     &       + ( 0.266)*r(527)  ! ISOP + O3 =
     &       + ( 0.728)*r(531)  ! APIN + O3 =
     &       + ( 0.500)*r(535)  ! ACYE + O3 =
     &       + ( 0.128)*r(549)  ! OLE1 + O3 =
     &       + ( 0.443)*r(553)  ! OLE2 + O3 =
     &       + ( 0.499)*r(559)  ! TERP + O3 =
     &       + ( 0.499)*r(563)  ! SESQ + O3 =
      newOH_O3 = PA(nn)
c
c --- OH directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newOH_phot'
      cpadesc(nn) = 'OH production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 23)  ! HONO =
     &       +          r( 28)  ! HNO3 =
     &       + ( 2.000)*r( 41)  ! H2O2 =
     &       +          r(223)  ! COOH =
     &       +          r(225)  ! ROOH =
     &       +          r(227)  ! R6PX =
     &       +          r(229)  ! RAPX =
     &       + ( 0.330)*r(258)  ! MACR =
     &       + ( 0.178)*r(278)  ! ACRO =
     &       +          r(280)  ! CO3H =
     &       +          r(282)  ! RO3H =
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
      PA(nn) =          r( 23)  ! HONO =
c
c --- OH loss
c
      nn = nn + 1
      ptname(nn)  = 'OH_loss'
      cpadesc(nn) = 'OH loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 22)  ! OH + NO =
     &       +          r( 24)  ! OH + HONO =
     &       +          r( 25)  ! OH + NO2 =
     &       +          r( 26)  ! OH + NO3 =
     &       +          r( 27)  ! OH + HNO3 =
     &       +          r( 29)  ! OH + CO =
     &       +          r( 30)  ! OH + O3 =
     &       +          r( 35)  ! PNA + OH =
     &       +          r( 42)  ! H2O2 + OH =
     &       +          r( 43)  ! OH + HO2 =
     &       +          r( 44)  ! OH + SO2 =
     &       +          r( 45)  ! OH + H2 =
     &       +          r(206)  ! HCHO + OH =
     &       +          r(208)  ! CCHO + OH =
     &       +          r(211)  ! RCHO + OH =
     &       +          r(214)  ! ACET + OH =
     &       +          r(216)  ! MEK + OH =
     &       +          r(218)  ! MEOH + OH =
     &       +          r(219)  ! FACD + OH =
     &       +          r(220)  ! AACD + OH =
     &       +          r(221)  ! PACD + OH =
     &       + ( 0.700)*r(222)  ! COOH + OH =
     &       + ( 0.256)*r(224)  ! ROOH + OH =
     &       + ( 0.160)*r(226)  ! R6PX + OH =
     &       + ( 0.861)*r(228)  ! RAPX + OH =
     &       +          r(232)  ! GLY + OH =
     &       +          r(235)  ! MGLY + OH =
     &       +          r(238)  ! CRES + OH =
     &       +          r(240)  ! NPHE + OH =
     &       +          r(243)  ! BALD + OH =
     &       +          r(246)  ! AFG1 + OH =
     &       +          r(249)  ! AFG2 + OH =
     &       +          r(252)  ! AFG3 + OH =
     &       +          r(254)  ! MACR + OH =
     &       +          r(259)  ! MVK + OH =
     &       +          r(263)  ! IPRD + OH =
     &       +          r(267)  ! PRD2 + OH =
     &       +          r(269)  ! RNO3 + OH =
     &       +          r(271)  ! GLYD + OH =
     &       +          r(274)  ! ACRO + OH =
     &       +          r(279)  ! CO3H + OH =
     &       +          r(281)  ! RO3H + OH =
     &       +          r(513)  ! CH4 + OH =
     &       +          r(514)  ! ETHE + OH =
     &       +          r(518)  ! PRPE + OH =
     &       +          r(522)  ! BD13 + OH =
     &       +          r(526)  ! ISOP + OH =
     &       +          r(530)  ! APIN + OH =
     &       + ( 0.300)*r(534)  ! ACYE + OH =
     &       + ( 0.884)*r(536)  ! BENZ + OH =
     &       + ( 0.688)*r(537)  ! TOLU + OH =
     &       + ( 0.761)*r(538)  ! MXYL + OH =
     &       + ( 0.802)*r(539)  ! OXYL + OH =
     &       + ( 0.722)*r(540)  ! PXYL + OH =
     &       + ( 0.770)*r(541)  ! B124 + OH =
     &       +          r(542)  ! ETOH + OH =
     &       +          r(543)  ! ALK1 + OH =
     &       +          r(544)  ! ALK2 + OH =
     &       +          r(545)  ! ALK3 + OH =
     &       +          r(546)  ! ALK4 + OH =
     &       +          r(547)  ! ALK5 + OH =
     &       +          r(548)  ! OLE1 + OH =
     &       +          r(552)  ! OLE2 + OH =
     &       + ( 0.798)*r(556)  ! ARO1 + OH =
     &       + ( 0.822)*r(557)  ! ARO2 + OH =
     &       +          r(558)  ! TERP + OH =
     &       +          r(562)  ! SESQ + OH =
      OH_loss = PA(nn)
c
c --- OH with CO
c
      nn = nn + 1
      ptname(nn)  = 'OHwCO'
      cpadesc(nn) = 'OH reaction rate with CO'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 29)  ! OH + CO =
c
c --- OH with isoprene
c
      nn = nn + 1
      ptname(nn)  = 'OHwISOP'
      cpadesc(nn) = 'OH reaction rate with isoprene'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(526)  ! ISOP + OH =
c
c --- OH with VOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwVOC'
      cpadesc(nn) = 'OH reaction rate with all VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(206)  ! HCHO + OH =
     &       +          r(208)  ! CCHO + OH =
     &       +          r(211)  ! RCHO + OH =
     &       +          r(214)  ! ACET + OH =
     &       +          r(216)  ! MEK + OH =
     &       +          r(218)  ! MEOH + OH =
     &       +          r(514)  ! ETHE + OH =
     &       +          r(518)  ! PRPE + OH =
     &       +          r(522)  ! BD13 + OH =
     &       +          r(526)  ! ISOP + OH =
     &       +          r(530)  ! APIN + OH =
     &       + ( 0.300)*r(534)  ! ACYE + OH =
     &       + ( 0.884)*r(536)  ! BENZ + OH =
     &       + ( 0.688)*r(537)  ! TOLU + OH =
     &       + ( 0.761)*r(538)  ! MXYL + OH =
     &       + ( 0.802)*r(539)  ! OXYL + OH =
     &       + ( 0.722)*r(540)  ! PXYL + OH =
     &       + ( 0.770)*r(541)  ! B124 + OH =
     &       +          r(542)  ! ETOH + OH =
     &       +          r(543)  ! ALK1 + OH =
     &       +          r(544)  ! ALK2 + OH =
     &       +          r(545)  ! ALK3 + OH =
     &       +          r(546)  ! ALK4 + OH =
     &       +          r(547)  ! ALK5 + OH =
     &       +          r(548)  ! OLE1 + OH =
     &       +          r(552)  ! OLE2 + OH =
     &       + ( 0.798)*r(556)  ! ARO1 + OH =
     &       + ( 0.822)*r(557)  ! ARO2 + OH =
     &       +          r(558)  ! TERP + OH =
     &       +          r(562)  ! SESQ + OH =
c
c --- OH with HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'OHwHRVOC'
      cpadesc(nn) = 'OH reaction rate with HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(514)  ! ETHE + OH =
     &       +          r(518)  ! PRPE + OH =
     &       +          r(522)  ! BD13 + OH =
     &       +          r(548)  ! OLE1 + OH =
     &       +          r(552)  ! OLE2 + OH =
c
c --- OH with Aromatics
c
      nn = nn + 1
      ptname(nn)  = 'OHwArom'
      cpadesc(nn) = 'OH reaction rate with aromatic VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.884)*r(536)  ! BENZ + OH =
     &       + ( 0.688)*r(537)  ! TOLU + OH =
     &       + ( 0.761)*r(538)  ! MXYL + OH =
     &       + ( 0.802)*r(539)  ! OXYL + OH =
     &       + ( 0.722)*r(540)  ! PXYL + OH =
     &       + ( 0.770)*r(541)  ! B124 + OH =
     &       + ( 0.798)*r(556)  ! ARO1 + OH =
     &       + ( 0.822)*r(557)  ! ARO2 + OH =
c
c --- OH with Alkanes (except methane)
c
      nn = nn + 1
      ptname(nn)  = 'OHwAlkane'
      cpadesc(nn) = 'OH reaction rate with alkanes (except methane)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(543)  ! ALK1 + OH =
     &       +          r(544)  ! ALK2 + OH =
     &       +          r(545)  ! ALK3 + OH =
     &       +          r(546)  ! ALK4 + OH =
     &       +          r(547)  ! ALK5 + OH =
c
c --- Total HCHO production
c
      nn = nn + 1
      ptname(nn)  = 'HCHO_prod'
      cpadesc(nn) = 'Total HCHO production rate'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 46)  ! MEO2 + NO =
     &       +          r( 48)  ! MEO2 + HO2 =
     &       +          r( 49)  ! MEO2 + NO3 =
     &       +          r( 50)  ! MEO2 + MEO2 =
     &       + ( 2.000)*r( 51)  ! MEO2 + MEO2 =
     &       + ( 0.750)*r( 55)  ! RO2C + MEO2 =
     &       + ( 0.750)*r( 60)  ! RO2X + MEO2 =
     &       +          r( 69)  ! MCO3 + MEO2 =
     &       +          r( 79)  ! RCO3 + MEO2 =
     &       +          r( 90)  ! BZC3 + MEO2 =
     &       + ( 0.400)*r( 98)  ! MPAN =
     &       +          r( 99)  ! MAC3 + NO =
     &       +          r(101)  ! MAC3 + NO3 =
     &       + ( 2.000)*r(102)  ! MAC3 + MEO2 =
     &       +          r(103)  ! MAC3 + RO2C =
     &       +          r(104)  ! MAC3 + RO2X =
     &       +          r(105)  ! MAC3 + MCO3 =
     &       +          r(106)  ! MAC3 + RCO3 =
     &       +          r(107)  ! MAC3 + BZC3 =
     &       + ( 2.000)*r(108)  ! MAC3 + MAC3 =
     &       +          r(218)  ! MEOH + OH =
     &       + ( 0.300)*r(222)  ! COOH + OH =
     &       +          r(223)  ! COOH =
     &       +          r(231)  ! GLY =
     &       + ( 0.100)*r(255)  ! MACR + O3 =
     &       + ( 0.340)*r(258)  ! MACR =
     &       + ( 0.050)*r(260)  ! MVK + O3 =
     &       + ( 0.124)*r(264)  ! IPRD + O3 =
     &       + ( 0.300)*r(266)  ! IPRD =
     &       + ( 0.002)*r(267)  ! PRD2 + OH =
     &       + ( 0.074)*r(270)  ! RNO3 =
     &       +          r(272)  ! GLYD =
     &       + ( 0.500)*r(275)  ! ACRO + O3 =
     &       + ( 0.340)*r(278)  ! ACRO =
     &       +          r(283)  ! XHCH + NO =
     &       +          r(285)  ! XHCH + NO3 =
     &       + ( 0.500)*r(286)  ! XHCH + MEO2 =
     &       + ( 0.500)*r(287)  ! XHCH + RO2C =
     &       + ( 0.500)*r(288)  ! XHCH + RO2X =
     &       +          r(289)  ! XHCH + MCO3 =
     &       +          r(290)  ! XHCH + RCO3 =
     &       +          r(291)  ! XHCH + BZC3 =
     &       +          r(292)  ! XHCH + MAC3 =
     &       +          r(515)  ! ETHE + O3 =
     &       + ( 0.500)*r(519)  ! PRPE + O3 =
     &       + ( 0.500)*r(523)  ! BD13 + O3 =
     &       + ( 0.400)*r(527)  ! ISOP + O3 =
     &       + ( 0.500)*r(549)  ! OLE1 + O3 =
     &       + ( 0.131)*r(553)  ! OLE2 + O3 =
     &       + ( 0.127)*r(559)  ! TERP + O3 =
     &       + ( 0.127)*r(563)  ! SESQ + O3 =
c
c --- HCHO production from HRVOC
c
      nn = nn + 1
      ptname(nn)  = 'nwHCHO_HRV'
      cpadesc(nn) = 'HCHO production rate from HRVOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(515)  ! ETHE + O3 =
     &       + ( 0.500)*r(519)  ! PRPE + O3 =
     &       + ( 0.500)*r(523)  ! BD13 + O3 =
     &       + ( 0.500)*r(549)  ! OLE1 + O3 =
     &       + ( 0.131)*r(553)  ! OLE2 + O3 =
c
c --- Total HO2 production, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'HO2_prod'
      cpadesc(nn) = 'Total HO2 production rate, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 26)  ! OH + NO3 =
     &       +          r( 29)  ! OH + CO =
     &       +          r( 30)  ! OH + O3 =
     &       +          r( 42)  ! H2O2 + OH =
     &       +          r( 44)  ! OH + SO2 =
     &       +          r( 45)  ! OH + H2 =
     &       +          r( 46)  ! MEO2 + NO =
     &       +          r( 49)  ! MEO2 + NO3 =
     &       + ( 2.000)*r( 51)  ! MEO2 + MEO2 =
     &       + ( 0.500)*r( 55)  ! RO2C + MEO2 =
     &       + ( 0.500)*r( 60)  ! RO2X + MEO2 =
     &       + ( 0.900)*r( 69)  ! MCO3 + MEO2 =
     &       +          r( 79)  ! RCO3 + MEO2 =
     &       +          r( 90)  ! BZC3 + MEO2 =
     &       +          r(102)  ! MAC3 + MEO2 =
     &       +          r(114)  ! XHO2 + NO =
     &       +          r(116)  ! XHO2 + NO3 =
     &       + ( 0.500)*r(117)  ! XHO2 + MEO2 =
     &       + ( 0.500)*r(118)  ! XHO2 + RO2C =
     &       + ( 0.500)*r(119)  ! XHO2 + RO2X =
     &       +          r(120)  ! XHO2 + MCO3 =
     &       +          r(121)  ! XHO2 + RCO3 =
     &       +          r(122)  ! XHO2 + BZC3 =
     &       +          r(123)  ! XHO2 + MAC3 =
     &       + ( 2.000)*r(204)  ! HCHO =
     &       +          r(206)  ! HCHO + OH =
     &       +          r(207)  ! HCHO + NO3 =
     &       +          r(209)  ! CCHO =
     &       +          r(212)  ! RCHO =
     &       +          r(218)  ! MEOH + OH =
     &       +          r(219)  ! FACD + OH =
     &       +          r(223)  ! COOH =
     &       +          r(225)  ! ROOH =
     &       + ( 0.142)*r(227)  ! R6PX =
     &       + ( 0.148)*r(228)  ! RAPX + OH =
     &       +          r(229)  ! RAPX =
     &       + ( 2.000)*r(230)  ! GLY =
     &       + ( 0.630)*r(232)  ! GLY + OH =
     &       + ( 0.630)*r(233)  ! GLY + NO3 =
     &       +          r(234)  ! MGLY =
     &       + ( 0.522)*r(247)  ! AFG1 + O3 =
     &       + ( 1.023)*r(248)  ! AFG1 =
     &       + ( 0.522)*r(250)  ! AFG2 + O3 =
     &       + ( 0.554)*r(253)  ! AFG3 + O3 =
     &       + ( 0.108)*r(255)  ! MACR + O3 =
     &       + ( 0.670)*r(258)  ! MACR =
     &       + ( 0.064)*r(260)  ! MVK + O3 =
     &       + ( 0.400)*r(264)  ! IPRD + O3 =
     &       + ( 1.233)*r(266)  ! IPRD =
     &       + ( 0.472)*r(267)  ! PRD2 + OH =
     &       + ( 0.189)*r(269)  ! RNO3 + OH =
     &       + ( 0.344)*r(270)  ! RNO3 =
     &       + ( 2.000)*r(272)  ! GLYD =
     &       + ( 0.830)*r(275)  ! ACRO + O3 =
     &       + ( 1.066)*r(278)  ! ACRO =
     &       +          r(485)  ! ZRN3 + NO3 =
     &       + ( 0.500)*r(486)  ! ZRN3 + MEO2 =
     &       + ( 0.500)*r(487)  ! ZRN3 + RO2C =
     &       + ( 0.500)*r(488)  ! ZRN3 + RO2X =
     &       +          r(489)  ! ZRN3 + MCO3 =
     &       +          r(490)  ! ZRN3 + RCO3 =
     &       +          r(491)  ! ZRN3 + BZC3 =
     &       +          r(492)  ! ZRN3 + MAC3 =
     &       + ( 0.160)*r(515)  ! ETHE + O3 =
     &       + ( 0.800)*r(517)  ! ETHE + O3P =
     &       + ( 0.165)*r(519)  ! PRPE + O3 =
     &       + ( 0.080)*r(523)  ! BD13 + O3 =
     &       + ( 0.250)*r(525)  ! BD13 + O3P =
     &       + ( 0.066)*r(527)  ! ISOP + O3 =
     &       + ( 0.009)*r(531)  ! APIN + O3 =
     &       + ( 0.300)*r(534)  ! ACYE + OH =
     &       + ( 1.500)*r(535)  ! ACYE + O3 =
     &       + ( 0.570)*r(536)  ! BENZ + OH =
     &       + ( 0.181)*r(537)  ! TOLU + OH =
     &       + ( 0.159)*r(538)  ! MXYL + OH =
     &       + ( 0.161)*r(539)  ! OXYL + OH =
     &       + ( 0.159)*r(540)  ! PXYL + OH =
     &       + ( 0.022)*r(541)  ! B124 + OH =
     &       + ( 0.950)*r(542)  ! ETOH + OH =
     &       + ( 0.095)*r(549)  ! OLE1 + O3 =
     &       + ( 0.094)*r(553)  ! OLE2 + O3 =
     &       + ( 0.123)*r(556)  ! ARO1 + OH =
     &       + ( 0.077)*r(557)  ! ARO2 + OH =
     &       + ( 0.078)*r(559)  ! TERP + O3 =
     &       + ( 0.078)*r(563)  ! SESQ + O3 =
c
c --- HO2 from O3 reactions with VOC
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_O3'
      cpadesc(nn) = 'HO2 production rate from O3 + VOC'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 0.522)*r(247)  ! AFG1 + O3 =
     &       + ( 0.522)*r(250)  ! AFG2 + O3 =
     &       + ( 0.554)*r(253)  ! AFG3 + O3 =
     &       + ( 0.108)*r(255)  ! MACR + O3 =
     &       + ( 0.064)*r(260)  ! MVK + O3 =
     &       + ( 0.400)*r(264)  ! IPRD + O3 =
     &       + ( 0.830)*r(275)  ! ACRO + O3 =
     &       + ( 0.160)*r(515)  ! ETHE + O3 =
     &       + ( 0.165)*r(519)  ! PRPE + O3 =
     &       + ( 0.080)*r(523)  ! BD13 + O3 =
     &       + ( 0.066)*r(527)  ! ISOP + O3 =
     &       + ( 0.009)*r(531)  ! APIN + O3 =
     &       + ( 1.500)*r(535)  ! ACYE + O3 =
     &       + ( 0.095)*r(549)  ! OLE1 + O3 =
     &       + ( 0.094)*r(553)  ! OLE2 + O3 =
     &       + ( 0.078)*r(559)  ! TERP + O3 =
     &       + ( 0.078)*r(563)  ! SESQ + O3 =
      newHO2_O3 = PA(nn)
c
c --- HO2 directly from photolysis, excluding from PNA
c
      nn = nn + 1
      ptname(nn)  = 'newHO2_pht'
      cpadesc(nn) = 'HO2 production rate directly from photolysis, excluding from PNA'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r(204)  ! HCHO =
     &       +          r(209)  ! CCHO =
     &       +          r(212)  ! RCHO =
     &       +          r(223)  ! COOH =
     &       +          r(225)  ! ROOH =
     &       + ( 0.142)*r(227)  ! R6PX =
     &       +          r(229)  ! RAPX =
     &       + ( 2.000)*r(230)  ! GLY =
     &       +          r(234)  ! MGLY =
     &       + ( 1.023)*r(248)  ! AFG1 =
     &       + ( 0.670)*r(258)  ! MACR =
     &       + ( 1.233)*r(266)  ! IPRD =
     &       + ( 0.344)*r(270)  ! RNO3 =
     &       + ( 2.000)*r(272)  ! GLYD =
     &       + ( 1.066)*r(278)  ! ACRO =
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
     &       + ( 2.000)*r(204)  ! HCHO =
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
     &       +          r(  6)  ! O3P + NO2 =
     &       +          r(  8)  ! O3 + NO2 =
     &       +          r( 12)  ! N2O5 =
     &       +          r( 27)  ! OH + HNO3 =
     &       + ( 0.390)*r( 34)  ! PNA =
     &       + ( 0.400)*r( 65)  ! PAN =
     &       + ( 0.400)*r( 75)  ! PAN2 =
     &       + ( 0.400)*r( 86)  ! PBZN =
     &       + ( 0.400)*r( 98)  ! MPAN =
c
c --- NO3 production from N2O5
c
      nn = nn + 1
      ptname(nn)  = 'N2O5toNO3'
      cpadesc(nn) = 'NO3 production from N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 12)  ! N2O5 =
c
c --- NO3 loss
c
      nn = nn + 1
      ptname(nn)  = 'NO3_loss'
      cpadesc(nn) = 'NO3 loss'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(  9)  ! NO + NO3 =
     &       +          r( 11)  ! NO2 + NO3 =
     &       +          r( 15)  ! NO2 + NO3 =
     &       +          r( 16)  ! NO3 =
     &       +          r( 17)  ! NO3 =
     &       +          r( 26)  ! OH + NO3 =
     &       +          r( 39)  ! NO3 + HO2 =
     &       + ( 2.000)*r( 40)  ! NO3 + NO3 =
     &       +          r( 49)  ! MEO2 + NO3 =
     &       +          r( 54)  ! RO2C + NO3 =
     &       +          r( 59)  ! RO2X + NO3 =
     &       +          r( 68)  ! MCO3 + NO3 =
     &       +          r( 78)  ! RCO3 + NO3 =
     &       +          r( 89)  ! BZC3 + NO3 =
     &       +          r(101)  ! MAC3 + NO3 =
     &       +          r(207)  ! HCHO + NO3 =
     &       +          r(210)  ! CCHO + NO3 =
     &       +          r(213)  ! RCHO + NO3 =
     &       +          r(233)  ! GLY + NO3 =
     &       +          r(236)  ! MGLY + NO3 =
     &       +          r(239)  ! CRES + NO3 =
     &       +          r(245)  ! BALD + NO3 =
     &       +          r(256)  ! MACR + NO3 =
     &       +          r(265)  ! IPRD + NO3 =
     &       +          r(273)  ! GLYD + NO3 =
     &       +          r(276)  ! ACRO + NO3 =
     &       +          r(516)  ! ETHE + NO3 =
     &       +          r(520)  ! PRPE + NO3 =
     &       +          r(524)  ! BD13 + NO3 =
     &       +          r(528)  ! ISOP + NO3 =
     &       +          r(532)  ! APIN + NO3 =
     &       +          r(550)  ! OLE1 + NO3 =
     &       +          r(554)  ! OLE2 + NO3 =
     &       +          r(560)  ! TERP + NO3 =
     &       +          r(564)  ! SESQ + NO3 =
c
c --- NO3 to N2O5
c
      nn = nn + 1
      ptname(nn)  = 'NO3toN2O5'
      cpadesc(nn) = 'NO3 conversion to N2O5'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 11)  ! NO2 + NO3 =
c
c --- Total production of ON compounds
c
      nn = nn + 1
      ptname(nn)  = 'ON_prod'
      cpadesc(nn) = 'Total production of organic nitrate (ON)'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 57)  ! RO2X + NO =
     &       +          r(135)  ! XNO2 + HO2 =
     &       + ( 0.500)*r(137)  ! XNO2 + MEO2 =
     &       + ( 0.500)*r(138)  ! XNO2 + RO2C =
     &       + ( 0.500)*r(139)  ! XNO2 + RO2X =
     &       + ( 0.500)*r(256)  ! MACR + NO3 =
     &       + ( 0.278)*r(265)  ! IPRD + NO3 =
     &       + ( 0.002)*r(276)  ! ACRO + NO3 =
     &       +          r(444)  ! XRN3 + HO2 =
     &       + ( 0.500)*r(446)  ! XRN3 + MEO2 =
     &       + ( 0.500)*r(447)  ! XRN3 + RO2C =
     &       + ( 0.500)*r(448)  ! XRN3 + RO2X =
     &       +          r(516)  ! ETHE + NO3 =
     &       +          r(520)  ! PRPE + NO3 =
     &       + ( 0.525)*r(524)  ! BD13 + NO3 =
     &       + ( 0.813)*r(528)  ! ISOP + NO3 =
     &       + ( 0.301)*r(532)  ! APIN + NO3 =
     &       + ( 0.226)*r(550)  ! OLE1 + NO3 =
     &       + ( 0.254)*r(554)  ! OLE2 + NO3 =
     &       + ( 0.485)*r(560)  ! TERP + NO3 =
     &       + ( 0.485)*r(564)  ! SESQ + NO3 =
     &       +          r(109)  ! TBUO + NO2 =
     &       +          r(443)  ! XRN3 + NO =
     &       +          r(445)  ! XRN3 + NO3 =
     &       + ( 0.500)*r(446)  ! XRN3 + MEO2 =
     &       + ( 0.500)*r(447)  ! XRN3 + RO2C =
     &       + ( 0.500)*r(448)  ! XRN3 + RO2X =
     &       +          r(449)  ! XRN3 + MCO3 =
     &       +          r(450)  ! XRN3 + RCO3 =
     &       +          r(451)  ! XRN3 + BZC3 =
     &       +          r(452)  ! XRN3 + MAC3 =
     &       +          r(483)  ! ZRN3 + NO =
     &       +          r(111)  ! BZO + NO2 =
      ON_prod = PA(nn)
c
c --- RNO3 production
c
      nn = nn + 1
      ptname(nn)  = 'RNO3_prod'
      cpadesc(nn) = 'Total production of the species RNO3'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(109)  ! TBUO + NO2 =
     &       +          r(443)  ! XRN3 + NO =
     &       +          r(445)  ! XRN3 + NO3 =
     &       + ( 0.500)*r(446)  ! XRN3 + MEO2 =
     &       + ( 0.500)*r(447)  ! XRN3 + RO2C =
     &       + ( 0.500)*r(448)  ! XRN3 + RO2X =
     &       +          r(449)  ! XRN3 + MCO3 =
     &       +          r(450)  ! XRN3 + RCO3 =
     &       +          r(451)  ! XRN3 + BZC3 =
     &       +          r(452)  ! XRN3 + MAC3 =
     &       +          r(483)  ! ZRN3 + NO =
c
c --- XN production
c
      nn = nn + 1
      ptname(nn)  = 'XN_prod'
      cpadesc(nn) = 'Total production of the species XN'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 57)  ! RO2X + NO =
     &       +          r(135)  ! XNO2 + HO2 =
     &       + ( 0.500)*r(137)  ! XNO2 + MEO2 =
     &       + ( 0.500)*r(138)  ! XNO2 + RO2C =
     &       + ( 0.500)*r(139)  ! XNO2 + RO2X =
     &       +          r(240)  ! NPHE + OH =
     &       +          r(242)  ! NPHE =
     &       + ( 0.500)*r(256)  ! MACR + NO3 =
     &       + ( 0.278)*r(265)  ! IPRD + NO3 =
     &       + ( 0.174)*r(269)  ! RNO3 + OH =
     &       + ( 0.002)*r(276)  ! ACRO + NO3 =
     &       +          r(444)  ! XRN3 + HO2 =
     &       + ( 0.500)*r(446)  ! XRN3 + MEO2 =
     &       + ( 0.500)*r(447)  ! XRN3 + RO2C =
     &       + ( 0.500)*r(448)  ! XRN3 + RO2X =
     &       +          r(516)  ! ETHE + NO3 =
     &       +          r(520)  ! PRPE + NO3 =
     &       + ( 0.525)*r(524)  ! BD13 + NO3 =
     &       + ( 0.813)*r(528)  ! ISOP + NO3 =
     &       + ( 0.301)*r(532)  ! APIN + NO3 =
     &       + ( 0.226)*r(550)  ! OLE1 + NO3 =
     &       + ( 0.254)*r(554)  ! OLE2 + NO3 =
     &       + ( 0.485)*r(560)  ! TERP + NO3 =
     &       + ( 0.485)*r(564)  ! SESQ + NO3 =
c
c --- HNO3 production
c
      nn = nn + 1
      ptname(nn)  = 'HNO3_prod'
      cpadesc(nn) = 'HNO3 production'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       + ( 2.000)*r( 13)  ! N2O5 + H2O =
     &       + ( 2.000)*r( 14)  ! N2O5 + H2O + H2O =
     &       +          r( 25)  ! OH + NO2 =
     &       + ( 0.200)*r( 39)  ! NO3 + HO2 =
     &       +          r(207)  ! HCHO + NO3 =
     &       +          r(210)  ! CCHO + NO3 =
     &       +          r(213)  ! RCHO + NO3 =
     &       +          r(233)  ! GLY + NO3 =
     &       +          r(236)  ! MGLY + NO3 =
     &       +          r(239)  ! CRES + NO3 =
     &       +          r(245)  ! BALD + NO3 =
     &       + ( 0.500)*r(256)  ! MACR + NO3 =
     &       + ( 0.150)*r(265)  ! IPRD + NO3 =
     &       +          r(273)  ! GLYD + NO3 =
     &       + ( 0.967)*r(276)  ! ACRO + NO3 =
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
     &       + ( 2.000)*r( 13)  ! N2O5 + H2O =
     &       + ( 2.000)*r( 14)  ! N2O5 + H2O + H2O =
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
     &       +          r( 28)  ! HNO3 =
     &       + ( 0.019)*r(269)  ! RNO3 + OH =
     &       +          r(270)  ! RNO3 =
c
c --- Net production of PAN compounds
c
      nn = nn + 1
      ptname(nn)  = 'PAN_prdNet'
      cpadesc(nn) = 'Net production of PAN species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r( 63)  ! MCO3 + NO2 =
     &       +          r( 73)  ! RCO3 + NO2 =
     &       +          r( 84)  ! BZC3 + NO2 =
     &       +          r( 96)  ! MAC3 + NO2 =
     &       -          r( 64)  ! PAN =
     &       -          r( 65)  ! PAN =
     &       -          r( 74)  ! PAN2 =
     &       -          r( 75)  ! PAN2 =
     &       -          r( 85)  ! PBZN =
     &       -          r( 86)  ! PBZN =
     &       -          r( 97)  ! MPAN =
     &       -          r( 98)  ! MPAN =
      PAN_prdNet = PA(nn)
c
c --- NO3 with VOC
c
      nn = nn + 1
      ptname(nn)  = 'NO3wVOC'
      cpadesc(nn) = 'NO3 reaction with VOC species'
      cpaunit(nn) = 'ppb hr-1'
      PA(nn) =
     &       +          r(207)  ! HCHO + NO3 =
     &       +          r(210)  ! CCHO + NO3 =
     &       +          r(213)  ! RCHO + NO3 =
     &       +          r(516)  ! ETHE + NO3 =
     &       +          r(520)  ! PRPE + NO3 =
     &       +          r(524)  ! BD13 + NO3 =
     &       +          r(528)  ! ISOP + NO3 =
     &       +          r(532)  ! APIN + NO3 =
     &       +          r(550)  ! OLE1 + NO3 =
     &       +          r(554)  ! OLE2 + NO3 =
     &       +          r(560)  ! TERP + NO3 =
     &       +          r(564)  ! SESQ + NO3 =
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

