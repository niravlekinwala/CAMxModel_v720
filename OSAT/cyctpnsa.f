c*** CYCTPNSA
c
      subroutine cyctpnsa(nspc,convfac,cold,delcon,rrxn_irr,
     &                    destRGN,prodRGN,NIT2RGN,RGN2NIT,
     &                    alphaRGN,betaRGN,gammaRGN,alphaNIT,cycNIT,
     &                    alphaNTR,betaNTR,alphaHN3,alphaTPN,betaTPN,
     &                    gammaTPN,cycTPN)
      use filunit
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine calculates the coefficients used in the chemistry
c     calculations of the nitrate species in SA.
c
c      Copyright 1996 - 2022
c     Ramboll 
c
c   Argument descriptions:
c     Inputs:
c       nspc       I  number of model species
c       convfac    R  conversion from ppm to umole/m3
c       cold       R  array of concentrations at last time step
c       delcon     R  array of change in concentrations
c       rrxn_irr   R  array of integrated reaction rates
c     Outputs:
c       destRGN    R  negative change in RGN (zero if RGN increasing)
c       prodRGN    R  positive change in RGN (zero if RGN decreasing)
c       NIT2RGN    R  conversion of NIT to RGN
c       RGN2NIT    R  conversion of RGN to NIT
c       alphaRGN   R  conversion of RGN to NTR
c       betaRGN    R  conversion of RGN to TPN (excluding cycTPN)
c       gammaRGN   R  conversion of RGN to NH3
c       alphaNIT   R  conversion of NIT to NTR
c       cycNIT     R  interconversion of NIT and RGN
c       alphaNTR   R  conversion of NTR to RGN
c       betaNTR    R  conversion of NTR to HN3
c       alphaHN3   R  conversion of HN3 to RGN
c       alphaTPN   R  conversion of TPN to HN3
c       betaTPN    R  conversion of TPN to RGN (excluding cycTPN)
c       gammaTPN   R  conversion of TPN to NTR
c       cycTPN     R  interconversion of TPN and RGN
c
c   Definitions of tracer families for reference:
c       NIT - sum of NO and HONO
c       RGN - sum of NO2, NO3, 2*N2O5 and INO3
c       TPN - sum of PAN compounds (and species that behave like PAN)
c       NIT - sum of organic nitrates (except those that behave like PAN)
c       HN3 - nitric acid (HNO3)
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     09/28/03   --gwilson--  Original development
c     12/29/06   --bkoo--     Revised for the updated SOA scheme
c     01/08/06   --bkoo--     Added Mechanism 6 (CB05)
c     05/17/07   --gwilson--  Changes to improve handling of small
c                             concentrations near lower bounds
c     12/22/13   --gwilson--  Added Mechanism 2 (CB6r2)
c     03/15/16   --gyarwood-  Changed for OSAT3 (NO2 apportioned)
c     06/24/16   --bkoo--     Updated Mechanism 4 (CB6r4)
c     10/16/17   --bkoo--     Changed DELCON size
c     11/27/17   --bkoo--     Updated for CB6r2h
c     01/16/18   --gyarwood-  Added SAPRC07T
c     02/21/20   --rlb--      Added Mechanism 1 (CB6r5)
c     03/24/21   --gy--       Added Mechanism 7 (CB7)
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer nspc
      real convfac
      real cold(nspc)
      real delcon(7,MXSPEC+1)
      real rrxn_irr(*)
      real destRGN
      real prodRGN
      real NIT2RGN
      real RGN2NIT
      real alphaRGN
      real betaRGN
      real gammaRGN
      real alphaNIT
      real cycNIT
      real alphaNTR
      real betaNTR
      real alphaHN3
      real alphaTPN
      real betaTPN
      real gammaTPN
      real cycTPN
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      real FUZZ
      parameter (FUZZ = 10.0)

      real lossPAN,lossPNA,lossPANX,lossPAN2,lossPBZN,lossMPAN,
     &     lossOPAN, lossINTR
      real prodPAN,prodPNA,prodPANX,prodPAN2,prodPBZN,prodMPAN,
     &     prodOPAN, prodINTR
      real cycPAN, cycPNA, cycPANX, cycPAN2, cycPBZN, cycMPAN,
     &     cycOPAN, cycINTR
      real sumRGN, sumTPN
      real delRGN, delNIT, delTPN
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      alphaNTR = 0.0
      betaNTR = 0.0
      alphaTPN = 0.0
      gammaTPN = 0.0
      betaTPN = 0.0
      cycTPN = 0.0
c
c  --- CB6r2h
c
      if( idmech .EQ. 3 ) then

        sumRGN = cold(kno2) + cold(kno3) + 2.*cold(kn2o5) + cold(kino3)
        delRGN = delcon(2,kno2) + delcon(2,kno3) + 2.*delcon(2,kn2o5)
     &                                                + delcon(2,kino3)
        delNIT = delcon(2,kno) + delcon(2,khono)
        cycNIT = MIN( 0.9999, 1. - EXP(-rrxn_irr(1)/sumRGN) )
        alphaNIT = MAX( 0.0, rrxn_irr(83) )

        if( (cold(kntr1) + cold(kntr2) + cold(kcron)) .GT. FUZZ*convfac*
     &                  (bdnl(kntr1) + bdnl(kntr2) + bdnl(kcron)) ) then
           alphaNTR = rrxn_irr(92)
           betaNTR  = rrxn_irr(215)
        endif

        alphaHN3 = rrxn_irr(46) + rrxn_irr(47)

        if( cold(khno3) .LT. FUZZ*convfac*bdnl(khno3) ) alphaHN3 = 0.0

        gammaRGN = MAX( 0.0, delcon(4,kHNO3) - alphaTPN - betaNTR
     &                                                    + alphaHN3 )
        sumTPN = cold(kpan) + cold(kpanx)
     &           + cold(kopan) + cold(kpna) + cold(kintr)
        delTPN = delcon(2,kpan) + delcon(2,kpanx)
     &           + delcon(2,kopan) + delcon(2,kpna) + delcon(2,kintr)

        if( sumTPN .GT. FUZZ*convfac*(bdnl(kpan) + bdnl(kpanx)
     &                  + bdnl(kopan) + bdnl(kpna) + bdnl(kintr)) ) then
           alphaTPN = 0.0
           gammaTPN = 0.5 * rrxn_irr(213)
           lossPAN = cold(kpan) *
     &         (1. - EXP(-(rrxn_irr(55)+rrxn_irr(56))/cold(kpan)) )
           prodPAN = cold(kno2) * (1. - EXP(-rrxn_irr(54)/cold(kno2)) )
           cycPAN = 0.5 * MIN(lossPAN, prodPAN)

           lossPANX = cold(kpanx) *
     &         (1. - EXP(-(rrxn_irr(63)+rrxn_irr(64)+rrxn_irr(214))
     &                                                   /cold(kpanx)) )
           prodPANX = cold(kno2) * (1. - EXP(-rrxn_irr(62)/cold(kno2)) )
           cycPANX = 0.5 * MIN(lossPANX, prodPANX)

           lossOPAN = cold(kopan) *
     &         (1. - EXP(-(rrxn_irr(209)+rrxn_irr(213))/cold(kopan)) )
           prodOPAN = cold(kno2) * (1. - EXP(-rrxn_irr(208)/cold(kno2)) )
           cycOPAN = 0.5 * MIN(lossOPAN, prodOPAN)

           lossPNA = cold(kpna) *
     &         (1. - EXP(-(rrxn_irr(49)+rrxn_irr(50)+rrxn_irr(51))
     &                                                   /cold(kpna)) )
           prodPNA = cold(kno2) * (1. - EXP(-rrxn_irr(48)/cold(kno2)) )
           cycPNA = 0.5 * MIN(lossPNA, prodPNA)

           lossINTR = cold(kintr) *
     &         (1. - EXP(-0.894*rrxn_irr(170)/cold(kintr)) )
           prodINTR = cold(kno2) * (1. - EXP(-0.1*rrxn_irr(151)/cold(kno2)) )
           cycINTR = 0.5 * AMIN1(lossINTR, prodINTR)

           cycTPN = cycPAN * cold(kpan)/sumTPN +
     &                cycPANX * cold(kpanx)/sumTPN +
     &                  cycOPAN * cold(kopan)/sumTPN +
     &                    cycPNA * cold(kpna)/sumTPN +
     &                      cycINTR * cold(kintr)/sumTPN

        endif

        alphaRGN = MAX( 0.0, delcon(2,kntr1) + delcon(2,kntr2) + delcon(2,kCRON)
     &                                       - alphaNIT - alphaTPN + betaNTR ) 
c
c  --- CB6r4 and CB6r5
c
      elseif( idmech .EQ. 1 .OR. idmech .EQ. 4 ) then

        sumRGN = cold(kno2) + cold(kno3) + 2.*cold(kn2o5) + cold(kino3)
        delRGN = delcon(2,kno2) + delcon(2,kno3) + 2.*delcon(2,kn2o5)
     &                                                + delcon(2,kino3)
        delNIT = delcon(2,kno) + delcon(2,khono)
        cycNIT = MIN( 0.9999, 1. - EXP(-rrxn_irr(1)/sumRGN) )
        alphaNIT = MAX( 0.0, rrxn_irr(83) )

        if( (cold(kntr1) + cold(kntr2) + cold(kcron)) .GT. FUZZ*convfac*
     &                  (bdnl(kntr1) + bdnl(kntr2) + bdnl(kcron)) ) then
           alphaNTR = rrxn_irr(92)
           betaNTR  = rrxn_irr(207)
        endif

        alphaHN3 = rrxn_irr(46) + rrxn_irr(47)

        if( cold(khno3) .LT. FUZZ*convfac*bdnl(khno3) ) alphaHN3 = 0.0

        gammaRGN = MAX( 0.0, delcon(4,kHNO3) - alphaTPN - betaNTR
     &                                                    + alphaHN3 )
        sumTPN = cold(kpan) + cold(kpanx)
     &           + cold(kopan) + cold(kpna) + cold(kintr)
        delTPN = delcon(2,kpan) + delcon(2,kpanx)
     &           + delcon(2,kopan) + delcon(2,kpna) + delcon(2,kintr)

        if( sumTPN .GT. FUZZ*convfac*(bdnl(kpan) + bdnl(kpanx)
     &                  + bdnl(kopan) + bdnl(kpna) + bdnl(kintr)) ) then
           alphaTPN = rrxn_irr(229)
           gammaTPN = 0.5 * rrxn_irr(205)
           lossPAN = cold(kpan) *
     &         (1. - EXP(-(rrxn_irr(55)+rrxn_irr(56))/cold(kpan)) )
           prodPAN = cold(kno2) * (1. - EXP(-rrxn_irr(54)/cold(kno2)) )
           cycPAN = 0.5 * MIN(lossPAN, prodPAN)

           lossPANX = cold(kpanx) *
     &         (1. - EXP(-(rrxn_irr(63)+rrxn_irr(64)+rrxn_irr(206))
     &                                                   /cold(kpanx)) )
           prodPANX = cold(kno2) * (1. - EXP(-rrxn_irr(62)/cold(kno2)) )
           cycPANX = 0.5 * MIN(lossPANX, prodPANX)

           lossOPAN = cold(kopan) *
     &         (1. - EXP(-(rrxn_irr(201)+rrxn_irr(205))/cold(kopan)) )
           prodOPAN = cold(kno2) * (1. - EXP(-rrxn_irr(200)/cold(kno2)) )
           cycOPAN = 0.5 * MIN(lossOPAN, prodOPAN)

           lossPNA = cold(kpna) *
     &         (1. - EXP(-(rrxn_irr(49)+rrxn_irr(50)+rrxn_irr(51))
     &                                                   /cold(kpna)) )
           prodPNA = cold(kno2) * (1. - EXP(-rrxn_irr(48)/cold(kno2)) )
           cycPNA = 0.5 * MIN(lossPNA, prodPNA)

           lossINTR = cold(kintr) *
     &         (1. - EXP(-0.894*rrxn_irr(163)/cold(kintr)) )
           prodINTR = cold(kno2) * (1. - EXP(-0.1*rrxn_irr(144)/cold(kno2)) )
           cycINTR = 0.5 * AMIN1(lossINTR, prodINTR)

           cycTPN = cycPAN * cold(kpan)/sumTPN +
     &                cycPANX * cold(kpanx)/sumTPN +
     &                  cycOPAN * cold(kopan)/sumTPN +
     &                    cycPNA * cold(kpna)/sumTPN +
     &                      cycINTR * cold(kintr)/sumTPN

        endif

        alphaRGN = MAX( 0.0, delcon(2,kntr1) + delcon(2,kntr2) + delcon(2,kCRON)
     &                                       - alphaNIT - alphaTPN + betaNTR ) 
c
c  --- CB7
c
      elseif( idmech .EQ. 7 ) then

        sumRGN = cold(kno2) + cold(kno3) + 2.*cold(kn2o5) + cold(kino3)
        delRGN = delcon(2,kno2) + delcon(2,kno3) + 2.*delcon(2,kn2o5)
     &                                                + delcon(2,kino3)
        delNIT = delcon(2,kno) + delcon(2,khono)
        cycNIT = MIN( 0.9999, 1. - EXP(-rrxn_irr(1)/sumRGN) )
        alphaNIT = MAX( 0.0, rrxn_irr(90) + 0.23*rrxn_irr(196) +
     &                      0.26*rrxn_irr(202) + 0.4*rrxn_irr(208) )

        if( (cold(kntr1) + cold(kntr2) + cold(kcron)) .GT. FUZZ*convfac*
     &                  (bdnl(kntr1) + bdnl(kntr2) + bdnl(kcron)) ) then
           alphaNTR = rrxn_irr(96) + rrxn_irr(174)
           betaNTR  = rrxn_irr(97)
        endif

        alphaHN3 = rrxn_irr(43) + rrxn_irr(44)

        if( cold(khno3) .LT. FUZZ*convfac*bdnl(khno3) ) alphaHN3 = 0.0

        gammaRGN = MAX( 0.0, delcon(4,kHNO3) - alphaTPN - betaNTR
     &                                                    + alphaHN3 )
        sumTPN = cold(kpan) + cold(kpanx)
     &           + cold(kopan) + cold(kpna) + cold(kintr)
        delTPN = delcon(2,kpan) + delcon(2,kpanx)
     &           + delcon(2,kopan) + delcon(2,kpna) + delcon(2,kintr)

        if( sumTPN .GT. FUZZ*convfac*(bdnl(kpan) + bdnl(kpanx)
     &                  + bdnl(kopan) + bdnl(kpna) + bdnl(kintr)) ) then
           alphaTPN = 0.0
           gammaTPN = 0.5 * rrxn_irr(72)
           lossPAN = cold(kpan) *
     &         (1. - EXP(-(rrxn_irr(58)+rrxn_irr(59))/cold(kpan)) )
           prodPAN = cold(kno2) * (1. - EXP(-rrxn_irr(57)/cold(kno2)) )
           cycPAN = 0.5 * MIN(lossPAN, prodPAN)

           lossPANX = cold(kpanx) *
     &         (1. - EXP(-(rrxn_irr(65)+rrxn_irr(66))/cold(kpanx)) )
           prodPANX = cold(kno2) * (1. - EXP(-rrxn_irr(64)/cold(kno2)) )
           cycPANX = 0.5 * MIN(lossPANX, prodPANX)

           lossOPAN = cold(kopan) *
     &         (1. - EXP(-(rrxn_irr(71)+0.5*rrxn_irr(72))/cold(kopan)) )
           prodOPAN = cold(kno2) * (1. - EXP(-rrxn_irr(70)/cold(kno2)) )
           cycOPAN = 0.5 * MIN(lossOPAN, prodOPAN)

           lossPNA = cold(kpna) *
     &         (1. - EXP(-(rrxn_irr(46)+rrxn_irr(47)+rrxn_irr(48))
     &                                                   /cold(kpna)) )
           prodPNA = cold(kno2) * (1. - EXP(-rrxn_irr(45)/cold(kno2)) )
           cycPNA = 0.5 * MIN(lossPNA, prodPNA)

           lossINTR = cold(kintr) *
     &         (1. - EXP(-0.5*rrxn_irr(194)/cold(kintr)) )
           prodINTR = cold(kno2) * (1. - EXP(-0.1*rrxn_irr(178)/cold(kno2)) )
           cycINTR = 0.5 * AMIN1(lossINTR, prodINTR)

           cycTPN = cycPAN * cold(kpan)/sumTPN +
     &                cycPANX * cold(kpanx)/sumTPN +
     &                  cycOPAN * cold(kopan)/sumTPN +
     &                    cycPNA * cold(kpna)/sumTPN +
     &                      cycINTR * cold(kintr)/sumTPN

        endif

        alphaRGN = MAX( 0.0, delcon(2,kntr1) + delcon(2,kntr2) + delcon(2,kCRON)
     &                                       - alphaNIT - alphaTPN + betaNTR ) 
c
c  --- SAPRC07T
c
      elseif( idmech .EQ. 5 ) then

        sumRGN = cold(kno2) + cold(kno3) + 2.*cold(kn2o5)
        delRGN = delcon(2,kno2) + delcon(2,kno3) + 2.*delcon(2,kn2o5)
        delNIT = delcon(2,kno) + delcon(2,khono)
        cycNIT = MIN( 0.9999, 1. - EXP(-rrxn_irr(1)/sumRGN) )
        alphaNIT = MAX( 0.0, (
     &              +      rrxn_irr(109)+      rrxn_irr(443)                     ! RNO3
     &              +      rrxn_irr(445)+0.500*rrxn_irr(446)+0.500*rrxn_irr(447)
     &              +0.500*rrxn_irr(448)+      rrxn_irr(449)+      rrxn_irr(450)
     &              +      rrxn_irr(451)+      rrxn_irr(452)
     &              +      rrxn_irr(111)                                         ! NPHE
     &              +      rrxn_irr( 57)+      rrxn_irr(135)+0.500*rrxn_irr(137) ! XN
     &              +0.500*rrxn_irr(138)+0.500*rrxn_irr(139)+      rrxn_irr(240)
     &              +      rrxn_irr(242)+0.500*rrxn_irr(256)+0.278*rrxn_irr(265)
     &              +      rrxn_irr(444)
     &              +0.500*rrxn_irr(446)+0.500*rrxn_irr(447)+0.500*rrxn_irr(448)
     &              +      rrxn_irr(516)+      rrxn_irr(520)
     &              +0.525*rrxn_irr(524)+0.813*rrxn_irr(528)+0.301*rrxn_irr(532)
     &              +0.226*rrxn_irr(550)+0.254*rrxn_irr(554)+0.485*rrxn_irr(560)
     &              +0.485*rrxn_irr(564)
     &                        ) )

        if( (cold(krno3) + cold(kxn)) .GT. FUZZ*convfac*
     &                  (bdnl(krno3) + bdnl(kxn)) ) then
           alphaNTR = 0.019*rrxn_irr(269)+rrxn_irr(270)
           betaNTR  = 0.0
        endif

        alphaHN3 = rrxn_irr(28) + rrxn_irr(27)

        if( cold(khno3) .LT. FUZZ*convfac*bdnl(khno3) ) alphaHN3 = 0.0

        gammaRGN = MAX( 0.0, delcon(4,kHNO3) - alphaTPN - betaNTR
     &                                                    + alphaHN3 )

        sumTPN = cold(kpan) + cold(kpan2)
     &           + cold(kmpan) + cold(kpbzn) + cold(kpna)
        delTPN = delcon(2,kpan) + delcon(2,kpan2)
     &           + delcon(2,kmpan) + delcon(2,kpbzn) + delcon(2,kpna)

        if( sumTPN .GT. FUZZ*convfac*(bdnl(kpan) + bdnl(kpan2)
     &                  + bdnl(kmpan) + bdnl(kpbzn) + bdnl(kpna)) ) then
           alphaTPN = 0.0
           gammaTPN = 0.0
           lossPAN = cold(kpan) *
     &         (1. - EXP(-(rrxn_irr(64)+rrxn_irr(65))/cold(kpan)) )
           prodPAN = cold(kno2) * (1. - EXP(-rrxn_irr(63)/cold(kno2)) )
           cycPAN = 0.5 * MIN(lossPAN, prodPAN)

           lossPAN2 = cold(kpan2) *
     &         (1. - EXP(-(rrxn_irr(74)+rrxn_irr(75)+rrxn_irr(214))
     &                                                   /cold(kpan2)) )
           prodPAN2 = cold(kno2) * (1. - EXP(-rrxn_irr(72)/cold(kno2)) )
           cycPAN2 = 0.5 * MIN(lossPAN2, prodPAN2)

           lossMPAN = cold(kmpan) *
     &         (1. - EXP(-(rrxn_irr(97)+rrxn_irr(98))/cold(kmpan)) )
           prodMPAN = cold(kno2) * (1. - EXP(-rrxn_irr(96)/cold(kno2)) )
           cycMPAN = 0.5 * MIN(lossMPAN, prodMPAN)

           lossPBZN = cold(kpbzn) *
     &         (1. - EXP(-(rrxn_irr(85)+rrxn_irr(86))/cold(kpbzn)) )
           prodPBZN = cold(kno2) * (1. - EXP(-rrxn_irr(84)/cold(kno2)) )
           cycPBZN = 0.5 * MIN(lossPBZN, prodPBZN)

           lossPNA = cold(kpna) *
     &         (1. - EXP(-(rrxn_irr(33)+rrxn_irr(34)+rrxn_irr(35))
     &                                                   /cold(kpna)) )
           prodPNA = cold(kno2) * (1. - EXP(-rrxn_irr(32)/cold(kno2)) )
           cycPNA = 0.5 * AMIN1(lossPNA, prodPNA)

           cycTPN = cycPAN * cold(kpan)/sumTPN +
     &                cycPAN2 * cold(kpan2)/sumTPN +
     &                  cycMPAN * cold(kmpan)/sumTPN +
     &                    cycPBZN * cold(kpbzn)/sumTPN +
     &                      cycPNA * cold(kpna)/sumTPN

        endif

        alphaRGN = MAX( 0.0, delcon(2,krno3) + delcon(2,kxn) + delcon(2,knphe)
     &                                       - alphaNIT - alphaTPN + betaNTR ) 
c
c  --- CB05 
c
      elseif( idmech .EQ. 6 ) then

        sumRGN = cold(kno2) + cold(kno3) + 2.*cold(kn2o5)
        delRGN = delcon(2,kno2) + delcon(2,kno3) + 2.*delcon(2,kn2o5)
        delNIT = delcon(2,kno) + delcon(2,khono)
        cycNIT = MIN( 0.9999, 1. - EXP(-rrxn_irr(1)/sumRGN) )
        alphaNIT = MAX( 0.0, rrxn_irr(55) )

        if( cold(kntr) .GT. FUZZ*convfac*bdnl(kntr) ) then
           alphaNTR = rrxn_irr(62)
           betaNTR  = rrxn_irr(61)
        endif

        alphaHN3 = rrxn_irr(29) + rrxn_irr(52)

        if( cold(khno3) .LT. FUZZ*convfac*bdnl(khno3) ) alphaHN3 = 0.0

        gammaRGN = MAX( 0.0, delcon(4,kHNO3) - alphaTPN - betaNTR
     &                                                    + alphaHN3 )
        sumTPN = cold(kpan) + cold(kpanx) + cold(kpna)
        delTPN = delcon(2,kpan) + delcon(2,kpanx) + delcon(2,kpna)

        if( sumTPN .GT. FUZZ*convfac*(bdnl(kpan) + bdnl(kpanx)
     &                  + bdnl(kpna)) ) then
           alphaTPN = 0.0
           gammaTPN = 0.0
           lossPAN = cold(kpan) *
     &         (1. - EXP(-(rrxn_irr(90)+rrxn_irr(91))/cold(kpan)) )
           prodPAN = cold(kno2) * (1. - EXP(-rrxn_irr(89)/cold(kno2)) )
           cycPAN = 0.5 * MIN(lossPAN, prodPAN)

           lossPANX = cold(kpanx) *
     &         (1. - EXP(-(rrxn_irr(105)+rrxn_irr(106)+rrxn_irr(107))
     &                                                    /cold(kpanx)) )
           prodPANX = cold(kno2) * (1. - EXP(-rrxn_irr(104)/cold(kno2)) )
           cycPANX = 0.5 * MIN(lossPANX, prodPANX)

           lossPNA = cold(kpna) *
     &         (1. - EXP(-(rrxn_irr(32)+rrxn_irr(33)+rrxn_irr(51))
     &                                                     /cold(kpna)) )
           prodPNA = cold(kno2) * (1. - EXP(-rrxn_irr(31)/cold(kno2)) )
           cycPNA = 0.5 * MIN(lossPNA, prodPNA)

           cycTPN = cycPAN * cold(kpan)/sumTPN +
     &                cycPANX * cold(kpanx)/sumTPN +
     &                  cycPNA * cold(kpna)/sumTPN
        endif

        alphaRGN = MAX( 0.0, delcon(2,kntr) 
     &                                       - alphaNIT - alphaTPN + betaNTR ) 
c
c  ---- stop if the mechanism is undefined -----
c
      else
         write(iout,'(//,a)') 'ERROR in CYCTPNSA:'
         write(iout,'(/,1X,A,I10)')
     &              'Unknown chemical mechanism ID number: ',idmech
         write(iout,'(/,1X,2A)')
     &              'Ozone source apportionment is not available ',
     &              'for this chemical mechanism'
         call camxerr()
      endif

c
c   --- set output variables that share definitions ---
c
      destRGN = MAX( 0., -delRGN )
      prodRGN = MAX( 0.,  delRGN )

      NIT2RGN = MAX( 0., -(delNIT + alphaNIT) )
      RGN2NIT = MAX( 0., delNIT + alphaNIT )

      betaTPN = MAX( 0., -delTPN )
      betaRGN = MAX( 0.,  delTPN )

      cycTPN  = MAX( 0., MIN( cycTPN,
     &                    0.33*(sumRGN-alphaRGN-betaRGN-gammaRGN),
     &                      0.33*(sumTPN-alphaTPN-betaTPN-gammaTPN) ) )

      betaTPN = betaTPN + cycTPN
      betaRGN = betaRGN + cycTPN
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end

