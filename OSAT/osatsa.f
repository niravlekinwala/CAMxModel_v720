c*** OSATSA
c
      subroutine osatsa(numcol,numrow,numlay,igrid,icell,jcell,kcell,
     &                              prdO3N,prdO3V,desto3,nspc,delcon,
     &                                   convfac,cold,rrxn_irr,dtime)
      use grid
      use chmstry
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine makes the "chemistry" adjustments to the tracer 
c     species.  The adjustments are based on the differences in 
c     concentrations of the regular model species before and after
c     the regular model chemistry.  This is essentially an adjustment
c     for production or decay.
c
c      Copyright 1996 - 2022
c     Ramboll 
c
c   Argument descriptions:
c     Inputs:
c       numcol    I  number of columns in this slice
c       numrow    I  number of rows in this slice
c       numlay    I  number of layers in this slice
c       igrid     I  grid number
c       icell     I  the X grid location of current cell
c       jcell     I  the X grid location of current cell
c       kcell     I  the vertical grid location of current layer
c       prdO3N    R  ozone production attributed to NOx
c       prdO3V    R  ozone production attributed to VOC
c       desto3    R  ozone destruction
c       nspc      I  number of model species
c       delcon    R  array of change in each species contrntrations
c       convfac   R  conversion factor used for lower bound value
c       cold      R  array of concentrations at last time step
c       rrxn_irr  R  array of reactions from last chemistry step
c       dtime     R  change in time for current time step 
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     05/22/97   --gwilson--  Now re-calculates the sum of VOC after
c                             adjustment using wtkoh.
c     09/20/03   --gwilson--  Changed the individual species args to
c                             a vector
c     07/20/04   --gwilson--  Changed for OSAT2 (ozone destruction)
c     07/27/15   --gyarwood-  Changed for OSAT3 (NO2 apportioned)
c     09/11/15   --bkoo--     Merged APCA
c     08/25/16   --gyarwood-  Updated for CB6r4
c     11/09/16   --cemery--   Added Baker APCA point source override option
c     10/16/17   --bkoo--     Changed DELCON size
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
      integer   numcol
      integer   numrow
      integer   numlay
      integer   igrid
      integer   icell
      integer   jcell
      integer   kcell
      real      prdO3N
      real      prdO3V
      real      desto3
      integer   nspc
      real      delcon(7,nspc+1)
      real      convfac
      real      cold(*)
      real      rrxn_irr(*)
      real      dtime
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer*8 idxcel, idx
      integer*8 idx_NIT, idx_HN3, idx_TPN, idx_NTR, idx_OO
      integer   i, i_N, iHN3, iTPN, iNTR, iOO, iVOC
      integer   istr1,istr2

      real      ptc1d(MXTRSP)
      real      sumNIT, sumHN3, sumTPN, sumNTR, sumRGN
      real      sumOO, sumO3, sumkOH, sumMIR, delVOC
      real      destRGN, prodRGN, NIT2RGN, RGN2NIT
      real      alphaRGN, betaRGN, gammaRGN, alphaNIT, cycNIT
      real      alphaNTR, betaNTR, alphaHN3
      real      alphaTPN, betaTPN, gammaTPN, cycTPN
      real      pairNOX, prodOO, destOO, OOtoO3, pairOO

      real      bioNIT, bioVOC, facbio, O3used
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  ipsa3d(igrid)-1 + DBLE(icell) + DBLE(numcol) * 
     &                DBLE(jcell-1) + DBLE(numcol) * DBLE(numrow) * 
     &                                                    DBLE(kcell-1)
      istr1 = 4
      istr2 = 6
      if (lapcapt) then
        istr1 = 7
        istr2 = 9
      endif

      do i = 1, ntotsp
         idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
         ptc1d(i) = ptconc(idx)
      enddo
c
c  --- get the sum of the tracer species for N ----
c
      sumNIT = 0.
      bioNIT = 0.
      do i=iptcls(idxipt(ITRNIT)),nptcls(idxipt(ITRNIT))
         sumNIT = sumNIT + ptc1d(i)
         if( ptname(i)(istr1:istr2) .EQ. '001' ) bioNIT = bioNIT + ptc1d(i)
      enddo
      sumNIT = MAX(BNDLPT, sumNIT)
      bioNIT = MAX(BNDLPT, bioNIT)
c
c  --- get the sum of the tracer species for HN3 ----
c
      sumHN3 = 0.
      do i=iptcls(idxipt(ITRHN3)),nptcls(idxipt(ITRHN3))
         sumHN3 = sumHN3 + ptc1d(i)
      enddo
      sumHN3 = MAX(BNDLPT, sumHN3)
c
c  --- get the sum of the tracer species for TPN ----
c
      sumTPN = 0.
      do i=iptcls(idxipt(ITRTPN)),nptcls(idxipt(ITRTPN))
         sumTPN = sumTPN + ptc1d(i)
      enddo
      sumTPN = MAX(BNDLPT, sumTPN)
c
c  --- get the sum of the tracer species for NTR ----
c
      sumNTR = 0.
      do i=iptcls(idxipt(ITRNTR)),nptcls(idxipt(ITRNTR))
         sumNTR = sumNTR + ptc1d(i)
      enddo
      sumNTR = MAX(BNDLPT, sumNTR)
c
c  --- get the sum of the tracer species for RGN ----
c
      sumRGN = 0.
      do i=iptcls(idxipt(ITRRGN)),nptcls(idxipt(ITRRGN))
         sumRGN = sumRGN + ptc1d(i)
      enddo
      sumRGN = MAX(BNDLPT, sumRGN)
c
c  --- call routine to calculate the cycle coefficients ---
c
      call cyctpnsa(nspec,convfac,cold,delcon,rrxn_irr,
     &              destRGN,prodRGN,NIT2RGN,RGN2NIT,
     &              alphaRGN,betaRGN,gammaRGN,alphaNIT,cycNIT,
     &              alphaNTR,betaNTR,alphaHN3,alphaTPN,betaTPN,
     &              gammaTPN,cycTPN)
c
c  --- adjust RGN, NIT, HN3, TPN and NTR tracers ----
c
      i_N  = iptcls(idxipt(ITRNIT)) - 1
      iHN3 = iptcls(idxipt(ITRHN3)) - 1
      iTPN = iptcls(idxipt(ITRTPN)) - 1
      iNTR = iptcls(idxipt(ITRNTR)) - 1
c
      do i=iptcls(idxipt(ITRRGN)),nptcls(idxipt(ITRRGN))
         i_N  = i_N  + 1
         iHN3 = iHN3 + 1
         iTPN = iTPN + 1
         iNTR = iNTR + 1
c
         pairNOX = ptc1d(i) + ptc1d(i_N)
c
         idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(i-1)
         ptconc(idx) =
     &      sumRGN * ( (1.0 - cycNIT) * ptc1d(i) / sumRGN +
     &                 (cycNIT * pairNOX / (sumNIT + sumRGN)) )
     &      + NIT2RGN * ptc1d(i_N) / sumNIT
     &          - (RGN2NIT + alphaRGN + betaRGN +
     &                           gammaRGN) * ptc1d(i) / sumRGN
     &              + alphaNTR * ptc1d(iNTR) / sumNTR
     &                   + alphaHN3 * ptc1d(iHN3) / sumHN3
     &                        + betaTPN * ptc1d(iTPN) / sumTPN
c
         idx_NIT = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(i_N-1)
         ptconc(idx_NIT) =
     &      sumNIT * ( (1.0-cycNIT) * ptc1d(i_N) / sumNIT +
     &                 (cycNIT * pairNOX / (sumNIT + sumRGN)) )
     &      + RGN2NIT * ptc1d(i) / sumRGN
     &         - (NIT2RGN + alphaNIT) * ptc1d(i_N) / sumNIT
c
         idx_HN3 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(iHN3-1)
         ptconc(idx_HN3) = ptconc(idx_HN3)
     &      + gammaRGN * ptc1d(i) / sumRGN
     &         + betaNTR * ptc1d(iNTR) / sumNTR
     &            + alphaTPN * ptc1d(iTPN) / sumTPN
     &               - alphaHN3 * ptc1d(iHN3) / sumHN3
c
         idx_TPN = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                        DBLE(numlay) * DBLE(iTPN-1)
         ptconc(idx_TPN) = ptconc(idx_TPN)
     &      + betaRGN * ptc1d(i) / sumRGN
     &         - (alphaTPN + betaTPN + gammaTPN) * ptc1d(iTPN) / sumTPN
c
         idx_NTR = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(iNTR-1)
         ptconc(idx_NTR) = ptconc(idx_NTR)
     &      + alphaNIT * ptc1d(i_N) / sumNIT
     &         + alphaRGN * ptc1d(i) / sumRGN
     &            + gammaTPN * ptc1d(iTPN) / sumTPN
     &               - (alphaNTR + betaNTR) * ptc1d(iNTR) / sumNTR
c
         ptconc(idx)     = MAX(BNDLPT,ptconc(idx))
         ptconc(idx_NIT) = MAX(BNDLPT,ptconc(idx_NIT))
         ptconc(idx_HN3) = MAX(BNDLPT,ptconc(idx_HN3))
         ptconc(idx_TPN) = MAX(BNDLPT,ptconc(idx_TPN))
         ptconc(idx_NTR) = MAX(BNDLPT,ptconc(idx_NTR))
      enddo
c
c   --- if no O3 tracer, return here ---
c
      if( .NOT. lozone ) return
c
c  --- get the sum of the tracer species for O3N/O3V ----
c
      sumO3 = 0.
      do i=iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3V))
         sumO3 = sumO3 + ptc1d(i)
      enddo
      sumO3 = MAX(BNDLPT, sumO3)
c
c  --- get the sum of the tracer species for OON/OOV ----
c
      sumOO = 0.
      do i=iptcls(idxipt(ITROON)),nptcls(idxipt(ITROOV))
         sumOO = sumOO + ptc1d(i)
      enddo
      sumOO = MAX(BNDLPT, sumOO)
c
c   --- set the local variable for VOC species concs ---
c
      delVOC = 0.
      do i=1,nspec
        if( lvocsp(i) ) delVOC = delVOC + delcon(1,i) * crbnum(i)
      enddo
c
c   --- loop over the VOC tracer species, calculate reactivity ---
c
      sumkOH = 0.
      sumMIR = 0.
      bioVOC = 0.
      do i=iptcls(idxipt(ITRVOC)),nptcls(idxipt(ITRVOC))
         sumkOH = sumkOH + ptc1d(i) * wtkoh(i)
         sumMIR = sumMIR + ptc1d(i) * wtmir(i)
         if( ptname(i)(istr1:istr2) .EQ. '001' )
     &                        bioVOC = bioVOC + ptc1d(i) * wtmir(i)
      enddo
      sumkOH = MAX(BNDLPT, sumkOH)
      sumMIR = MAX(BNDLPT, sumMIR)
      bioVOC = MAX(BNDLPT, bioVOC)
c
c   ---- decay VOC tracers ---
c
      do i=iptcls(idxipt(ITRVOC)),nptcls(idxipt(ITRVOC))
         idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(i-1)
         ptconc(idx) = ptconc(idx) + 
     &                 delVOC * ptc1d(i) * wtkoh(i) / sumkOH
         ptconc(idx) = MAX(BNDLPT,ptconc(idx))
      enddo
c
c  --- update OON/OOV and allocate any ozone from NO2 destruction ---
c
      prodOO = prodRGN
      destOO = destRGN

      OOtoO3 = MIN( destOO, prdO3N + prdO3V )
      if( prdO3V .GT. 0.) then
         prdO3V = prdO3V - OOtoO3
      else
         prdO3N = prdO3N - OOtoO3
      endif
c
      iOO = iptcls(idxipt(ITROON)) - 1
      do i=iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3V))
         iOO = iOO + 1
         idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(i-1)
         idx_OO = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iOO-1)
         pairOO = ptc1d(i) + ptc1d(iOO)
         ptconc(idx) =
     &      sumO3 * ( (1.0 - cycNIT) * ptc1d(i) / sumO3 +
     &                 (cycNIT * pairOO / (sumO3 + sumOO)) )
     &         + OOtoO3 * ptc1d(iOO) / sumOO
         ptconc(idx_OO) =
     &      sumOO * ( (1.0-cycNIT) * ptc1d(iOO) / sumOO +
     &                 (cycNIT * pairOO / (sumO3 + sumOO)) )
     &         - destOO * ptc1d(iOO) / sumOO
     &            + prodOO * ptc1d(i) / sumO3
         ptconc(idx_OO) = MAX(BNDLPT,ptconc(idx_OO))
      enddo
c
c  --- allocate any ozone destruction across all tracers ---
c
      if( desto3 .LT. 0. ) then
         do i=iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3V))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + desto3 * ptc1d(i) / sumO3
         enddo
      endif
c
c  --- APCA
c
      if ( lapca ) then
c
c  --- set the biogenics contribution factor ---
c
         facbio = MIN( bioNIT/sumNIT, bioVOC/sumMIR )
c
c  --- allocate VOC sensitive ozone production based on MIRs,
c      increase the O3V tracers except where limited by APCA ---
c
         if( prdO3V .GT. 0. ) then
            O3used = 0.
            iVOC = iptcls(idxipt(ITRVOC)) - 1
            do i=iptcls(idxipt(ITRO3V)),nptcls(idxipt(ITRO3V))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(i-1)
               iVOC = iVOC + 1
c
c  --- add O3 to biogenics group (group 1) based on it's
c      mimimum potential contribution ---
c
               if( ptname(i)(istr1:istr2) .EQ. '001' ) then
                  ptconc(idx) = ptconc(idx) + facbio * prdO3V *
     &                          ptc1d(iVOC) * wtmir(iVOC) / bioVOC
                  O3used = O3used + facbio * prdO3V *
     &                          ptc1d(iVOC) * wtmir(iVOC) / bioVOC
c
c  --- add O3 to other groups based on contribution of VOC conc ---
c
               else
                  ptconc(idx) = ptconc(idx) + prdO3V *
     &                          ptc1d(iVOC) * wtmir(iVOC) / sumMIR
                  O3used = O3used + prdO3V *
     &                          ptc1d(iVOC) * wtmir(iVOC) / sumMIR
               endif
            enddo
c
c  --- add any ozone not accounted for to O3N tracers, based on
c      contribution to anthropognic NOx ---
c
            if( sumNIT .GT. bioNIT .AND. prdO3V .GT. O3used ) then
               i_N = iptcls(idxipt(ITRNIT)) - 1
               do i=iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3N))
                  idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
                  i_N = i_N + 1
                  if( ptname(i)(istr1:istr2) .NE. '001' ) then
                     ptconc(idx) = ptconc(idx) + (prdO3V - O3used) *
     &                             ptc1d(i_N) / (sumNIT - bioNIT)
                  endif
               enddo
            endif
         endif
c
c  --- allocate NOx sensitive ozone production,
c      increase the O3N tracers except where limited by APCA ---
c
         if( prdO3N .GT. 0. ) then
            O3used = 0.
            i_N = iptcls(idxipt(ITRNIT)) - 1
            do i=iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3N))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
               i_N = i_N + 1
c
c  --- add O3 to biogenics group (group 1) based on it's
c      mimimum potential contribution ---
c
               if( ptname(i)(istr1:istr2) .EQ. '001' ) then
                  ptconc(idx) = ptconc(idx) + facbio * prdO3N *
     &                                        ptc1d(i_N) / bioNIT
                  O3used = O3used + facbio * prdO3N *
     &                                        ptc1d(i_N) / bioNIT
c
c  --- add O3 to other groups based on contribution of NOx conc ---
c
               else
                  ptconc(idx) = ptconc(idx) + prdO3N *
     &                                        ptc1d(i_N) / sumNIT
                  O3used = O3used + prdO3N * ptc1d(i_N) / sumNIT
               endif
            enddo
c
c  --- add any ozone not accounted for to O3V tracers, based on
c      contribution to anthropognic VOC ---
c
            if( sumMIR .GT. bioVOC .AND. prdO3N .GT. O3used ) then
               iVOC = iptcls(idxipt(ITRVOC)) - 1
               do i=iptcls(idxipt(ITRO3V)),nptcls(idxipt(ITRO3V))
                  idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
                  iVOC = iVOC + 1
                  if( ptname(i)(istr1:istr2) .NE. '001' ) then
                     ptconc(idx) = ptconc(idx) + (prdO3N - O3used) *
     &                  ptc1d(iVOC) * wtmir(iVOC) / (sumMIR - bioVOC)
                  endif
               enddo
            endif
         endif
c
c  --- OSAT
c
      else
c
c  --- allocate VOC sensitive ozone production based on MIRs,
c      increase the O3V tracers ---
c
         if( prdO3V .GT. 0. ) then
            iVOC = iptcls(idxipt(ITRVOC)) - 1
            do i=iptcls(idxipt(ITRO3V)),nptcls(idxipt(ITRO3V))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
               iVOC = iVOC + 1
               ptconc(idx) = ptconc(idx) + prdO3V *
     &                       ptc1d(iVOC) * wtmir(iVOC) / sumMIR
            enddo
         endif
c
c  --- allocate NOx sensitive ozone production,
c      increase the O3N tracers ---
c
         if( prdO3N .GT. 0. ) then
            i_N = iptcls(idxipt(ITRNIT)) - 1
            do i=iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3N))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
               i_N = i_N + 1
               ptconc(idx) = ptconc(idx) + prdO3N *
     &                                     ptc1d(i_N) / sumNIT
            enddo
         endif
c
      endif ! APCA?
c
c  --- lower bound of ozone tracers
c
      do i=iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3V))
         idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(i-1)
         ptconc(idx) = MAX(BNDLPT,ptconc(idx))
      enddo
c
c  --- decay the timing tracer species --- 
c
      if(ntrtim .GT. 0 ) then
         do i=ipttim+1,nsaspc,2
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) * EXP( -0.08333333 * dtime )
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
