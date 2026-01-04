c*** PSATSA
c
      subroutine psatsa(numcol,numrow,numlay,igrid,icell,jcell,kcell,
     &                                         nspc,delcon,nprec,yfac)
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
c     species for PSAT.  The adjustments are based on the differences in 
c     concentrations of the regular model species before and after
c     the regular model chemistry.  This is essentially an adjustment
c     for production or decay.
c
c      Copyright 1996 - 2022
c     Ramboll 
c
c   Argument descriptions:
c     Inputs:
c       numcol  I number of columns in this slice
c       numrow  I number of rows in this slice
c       numlay  I number of layers in this slice
c       igrid     I  grid number
c       icell     I  the X grid location of current cell
c       jcell     I  the X grid location of current cell
c       kcell     I  the vertical grid location of current layer
c       nspc      I  number of model species
c       delcon    R  array of change in concentrations total
c       nprec     I  number of SOA precursors
c       yfac      R  array of weighting factors for SOA
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     09/28/03   --gwilson--  Original development
c     12/29/06   --bkoo--     Revised for the updated SOA scheme
c     05/04/07   --gwilson--  Changed call to cyctpnsa so that conversion
c                             factor is set for lower bound values -
c                             the conversion is passed here from chemdriv
c     03/18/14   --bkoo--     Revised for benzene SOA
c     08/25/16   --bkoo--     Updated for new SOAP
c     10/16/17   --bkoo--     Fixed DELARO that was missing IVOA changes
c                             Updated for SOA photolysis
c     01/12/18   --bkoo--     Removed BNZA/TOLA/XYLA/ISP/TRP
c     01/09/19   --cemery--   Added DMS
c     06/18/19   --bkoo--     Updated delARO/delTRP/delSQT for SAPRC07T
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
      integer   nspc
      real      delcon(7,nspc+1)
      integer   nprec
      real      yfac(3,nprec)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   icg1, icg2, isopa, iaro, iisp, itrp, isqt
      integer   iso2, idms, ihg0, i
      integer*8 idxcel, idx, jdx, kdx, ldx, idx0, idx2
      real      delARO, delISP, delTRP, delSQT, delso2, delps4, deldms, delhg2
      real      delCG1, delCG2, delNV1
      real      sumaro, sumisp, sumtrp, sumsqt, sumso2, sumps4, sumdms
      real      sumhg0, sumhg2

      real, parameter :: wtsopb = 220.0 ! must be consistent with MWSOPB in SOAP.INC
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  DBLE(ipsa3d(igrid)-1) + DBLE(icell) + 
     &            DBLE(numcol) * DBLE(jcell-1) + 
     &                 DBLE(numcol) * DBLE(numrow) * DBLE(kcell-1)
c
c----------------------------------------------------------------------
c  ---  SOA tracers ----
c----------------------------------------------------------------------
c
      if( lsoa ) then
         delARO = delcon(1,kbenz)                   ! [umol/m3] for gases
     &          + delcon(1,ktol)  + delcon(1,ktolu) + delcon(1,karo1)
     &          + delcon(1,kxyl)  + delcon(1,koxyl) + delcon(1,kmxyl)
     &          + delcon(1,kpxyl) + delcon(1,kb124) + delcon(1,karo2)
     &          + delcon(1,kivoa)
         delISP = delcon(1,kisop)
         delTRP = delcon(1,kterp) + delcon(1,kapin)
         delSQT = delcon(1,ksqt)  + delcon(1,ksesq)
         delCG1 = delcon(1,kcg1)
         delCG2 = delcon(1,kcg2)
         delNV1 = delcon(3,ksopa) - delcon(5,ksopa) ! adjust non-vol SOA delta [ug/m3] to exclude
     &                            - delcon(6,ksopa) ! changes due to SOA polymerization/photolysis
c
c  --- get the sum of the tracer species for ARO/CG1 yield ---
c
         icg1 = iptcls(idxipt(ITRCG1)) - 1
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(i-1)
            icg1 = icg1 + 1
            sumaro = sumaro + ptconc(idx) *
     &                      ( yhrates(icg1)*yfac(1,1)
     &                      + ylrates(icg1)*yfac(2,1) )
         enddo
c
c  --- make adjustment to CG1 tracers based on change in CG1 and
c      yield rates for ARO ---
c
         iaro = iptcls(idxipt(ITRARO)) - 1
         do i=iptcls(idxipt(ITRCG1)),nptcls(idxipt(ITRCG1))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
            iaro = iaro + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(iaro-1)
            ptconc(idx) = ptconc(idx) + delCG1 * ( ptconc(jdx) *
     &                                  ( yhrates(i)*yfac(1,1)
     &                                  + ylrates(i)*yfac(2,1) ) )
     &                                                    / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ARO/CG2 yield ---
c
         icg2 = iptcls(idxipt(ITRCG2)) - 1
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
            icg2 = icg2 + 1
            sumaro = sumaro + ptconc(idx) *
     &                      ( yhrates(icg2)*yfac(1,1)
     &                      + ylrates(icg2)*yfac(2,1) )
         enddo
c
c  --- make adjustment to CG2 tracers based on change in CG2 and
c      yield rates for ARO ---
c
         iaro = iptcls(idxipt(ITRARO)) - 1
         do i=iptcls(idxipt(ITRCG2)),nptcls(idxipt(ITRCG2))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            iaro = iaro + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iaro-1)
            ptconc(idx) = ptconc(idx) + delCG2 * ( ptconc(jdx) *
     &                                  ( yhrates(i)*yfac(1,1)
     &                                  + ylrates(i)*yfac(2,1) ) )
     &                                                    / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ARO/SOPA yield ---
c      SPECIAL CASE: non-volatile CG; ARO -> SOPA directly (skip CG)
c
         isopa = iptcls(idxipt(ITRPPA)) - 1
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            isopa = isopa + 1
            sumaro = sumaro + ptconc(idx) *
     &                      ( yhrates(isopa)*yfac(1,1)
     &                      + ylrates(isopa)*yfac(2,1) )
         enddo
c
c  --- make adjustment to SOPA tracers based on change in SOPA and
c      yield rates for ARO ---
c      SPECIAL CASE: non-volatile CG; ARO -> SOPA directly (skip CG)
c
         iaro = iptcls(idxipt(ITRARO)) - 1
         do i=iptcls(idxipt(ITRPPA)),nptcls(idxipt(ITRPPA))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            iaro = iaro + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iaro-1)
            ptconc(idx) = ptconc(idx) + delNV1 * ( ptconc(jdx) *
     &                                  ( yhrates(i)*yfac(1,1)
     &                                  + ylrates(i)*yfac(2,1) ) )
     &                                                    / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ARO ---
c
         sumaro = 0.
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            sumaro = sumaro + ptconc(idx)*wtkoh(i)
         enddo
c
c  --- make adjustment to ARO based on change in aromatic precursors and kOH ---
c
         do i=iptcls(idxipt(ITRARO)),nptcls(idxipt(ITRARO))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delARO * 
     &                                  ptconc(idx)*wtkoh(i) / sumaro
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- get the sum of the tracer species for ISP/TRP/SQT ---
c
         sumisp = 0.
         do i=iptcls(idxipt(ITRISP)),nptcls(idxipt(ITRISP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            sumisp = sumisp + ptconc(idx)
         enddo

         sumtrp = 0.
         do i=iptcls(idxipt(ITRTRP)),nptcls(idxipt(ITRTRP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(i-1)
            sumtrp = sumtrp + ptconc(idx)
         enddo

         sumsqt = 0.
         do i=iptcls(idxipt(ITRSQT)),nptcls(idxipt(ITRSQT))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                              DBLE(numlay) * DBLE(i-1)
            sumsqt = sumsqt + ptconc(idx)
         enddo
c
c  --- make adjustment to CG3 tracers based on change in ISP/TRP/SQT and
c      distribution of ISP/TRP/SQT ----
c
         iisp = iptcls(idxipt(ITRISP)) - 1
         itrp = iptcls(idxipt(ITRTRP)) - 1
         isqt = iptcls(idxipt(ITRSQT)) - 1
         do i=iptcls(idxipt(ITRCG3)),nptcls(idxipt(ITRCG3))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(i-1)
            iisp = iisp + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(iisp-1)
            itrp = itrp + 1
            kdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(itrp-1)
            isqt = isqt + 1
            ldx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                             DBLE(numlay) * DBLE(isqt-1)
            ptconc(idx) = ptconc(idx) - delISP * yfac(1,5) *
     &                                           ptconc(jdx) / sumisp
     &                                - delTRP * yfac(1,6) *
     &                                           ptconc(kdx) / sumtrp
     &                                - delSQT * yfac(1,7) *
     &                                           ptconc(ldx) / sumsqt
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to CG4 tracers based on change in ISP/TRP/SQT and
c      distribution of ISP/TRP/SQT ----
c
         iisp = iptcls(idxipt(ITRISP)) - 1
         itrp = iptcls(idxipt(ITRTRP)) - 1
         isqt = iptcls(idxipt(ITRSQT)) - 1
         do i=iptcls(idxipt(ITRCG4)),nptcls(idxipt(ITRCG4))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(i-1)
            iisp = iisp + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(iisp-1)
            itrp = itrp + 1
            kdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(itrp-1)
            isqt = isqt + 1
            ldx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                            DBLE(numlay) * DBLE(isqt-1)
            ptconc(idx) = ptconc(idx) - delISP * yfac(2,5) *
     &                                           ptconc(jdx) / sumisp
     &                                - delTRP * yfac(2,6) *
     &                                           ptconc(kdx) / sumtrp
     &                                - delSQT * yfac(2,7) *
     &                                           ptconc(ldx) / sumsqt
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to SOPB tracers based on change in ISP/TRP/SQT and
c      distribution of ISP/TRP/SQT ----
c      SPECIAL CASE: non-volatile CG; ISP/TRP/SQT -> SOPB directly (skip CG)
c
         iisp = iptcls(idxipt(ITRISP)) - 1
         itrp = iptcls(idxipt(ITRTRP)) - 1
         isqt = iptcls(idxipt(ITRSQT)) - 1
         do i=iptcls(idxipt(ITRPPB)),nptcls(idxipt(ITRPPB))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            iisp = iisp + 1
            jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(iisp-1)
            itrp = itrp + 1
            kdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(itrp-1)
            isqt = isqt + 1
            ldx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(isqt-1)
            ptconc(idx) = ptconc(idx) - wtsopb * 
     &                                ( delISP * yfac(3,5) *
     &                                           ptconc(jdx) / sumisp
     &                                + delTRP * yfac(3,6) *
     &                                           ptconc(kdx) / sumtrp
     &                                + delSQT * yfac(3,7) *
     &                                           ptconc(ldx) / sumsqt )
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to ISP based on change in ISP ---
c
         do i=iptcls(idxipt(ITRISP)),nptcls(idxipt(ITRISP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                 DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delISP * ptconc(idx) / sumisp
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to TRP based on change in TRP ---
c
         do i=iptcls(idxipt(ITRTRP)),nptcls(idxipt(ITRTRP))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                     DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delTRP * ptconc(idx) / sumtrp
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to SQT based on change in SQT ---
c
         do i=iptcls(idxipt(ITRSQT)),nptcls(idxipt(ITRSQT))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delSQT * ptconc(idx) / sumsqt
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
      endif
c
c----------------------------------------------------------------------
c  ---  Sulfate tracers ---
c----------------------------------------------------------------------
c
      if( lsulfate ) then
         delso2 = delcon(3,kso2)
         delps4 = delcon(3,kpso4)
c
c  --- get the sum of the tracer species for SO2 ---
c
         sumso2 = 0.
         do i=iptcls(idxipt(ITRSO2)),nptcls(idxipt(ITRSO2))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumso2 = sumso2 + ptconc(idx)
         enddo
c
c  --- get the sum of the tracer species for DMS ---
c
         if( ldmschm ) then
            idms = iptcls(idxipt(ITRDMS)) - 1
            deldms = delcon(1,kdms)
            sumdms = 0.
            do i=iptcls(idxipt(ITRDMS)),nptcls(idxipt(ITRDMS))
               idx = idxcel + numcol * numrow * numlay * (i-1)
               sumdms = sumdms + ptconc(idx)
            enddo
         endif
c
c  --- if PSO4 change is positive, make adjustment to PS4 
c      tracers based on distribtion of SO2 ---
c
         if( delps4 .GT. 0. ) then
            iso2 = iptcls(idxipt(ITRSO2)) - 1
            do i=iptcls(idxipt(ITRPS4)),nptcls(idxipt(ITRPS4))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
               iso2 = iso2 + 1
               jdx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                               DBLE(numlay) * DBLE(iso2-1)
               ptconc(idx) = ptconc(idx) + delps4 * ptconc(jdx) / sumso2
               ptconc(idx) = MAX(BNDLPT,ptconc(idx))
            enddo
c
c  --- if PSO4 change is negative, make adjustment to PS4 
c      tracers based on distribtion of PS4 ---
c
         else
c
c  --- get the sum of the tracer species for PS4 ---
c
            sumps4 = 0.
            do i=iptcls(idxipt(ITRPS4)),nptcls(idxipt(ITRPS4))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
               sumps4 = sumps4 + ptconc(idx)
            enddo
            do i=iptcls(idxipt(ITRPS4)),nptcls(idxipt(ITRPS4))
               idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                    DBLE(numlay) * DBLE(i-1)
               ptconc(idx) = ptconc(idx) + delps4 * ptconc(idx) / sumps4
               ptconc(idx) = MAX(BNDLPT,ptconc(idx))
            enddo
         endif
c
c  --- make adjustment to SO2 tracers based on change in SO2 (and DMS) ---
c
         do i=iptcls(idxipt(ITRSO2)),nptcls(idxipt(ITRSO2))
            idx = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                  DBLE(numlay) * DBLE(i-1)
            ptconc(idx) = ptconc(idx) + delso2 * ptconc(idx) / sumso2
            if( ldmschm ) then
               idms = idms + 1
               jdx = idxcel + numcol * numrow * numlay * (idms-1)
               ptconc(idx) = ptconc(idx) - deldms * ptconc(jdx) / sumdms
            endif
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
         enddo
c
c  --- make adjustment to DMS tracers based on change in DMS ---
c
         if( ldmschm ) then
           do i=iptcls(idxipt(ITRDMS)),nptcls(idxipt(ITRDMS))
              idx = idxcel + numcol * numrow * numlay * (i-1)
              ptconc(idx) = ptconc(idx) + deldms * ptconc(idx) / sumdms
              ptconc(idx) = MAX(BNDLPT,ptconc(idx))
           enddo
         endif
      endif
c
c----------------------------------------------------------------------
c  ---  Mercury tracers ----
c----------------------------------------------------------------------
c
      if( lmercury ) then
         delhg2 = delcon(3,khg2)

c  --- get the sum of the tracer species ---
c
         sumhg0 = 0.
         do i=iptcls(idxipt(ITRHG0)),nptcls(idxipt(ITRHG0))
            idx0 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumhg0 = sumhg0 + ptconc(idx0)
         enddo
         sumhg2 = 0.
         do i=iptcls(idxipt(ITRHG2)),nptcls(idxipt(ITRHG2))
            idx2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
            sumhg2 = sumhg2 + ptconc(idx2)
         enddo
c
c  --- HG0 being oxidized to HG2 ----
c
         if( delhg2 .GT. 0 ) then
            ihg0 = iptcls(idxipt(ITRHG0)) - 1
            do i=iptcls(idxipt(ITRHG2)),nptcls(idxipt(ITRHG2))
               idx2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(i-1)
               ihg0 = ihg0 + 1
               idx0 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                DBLE(numlay) * DBLE(ihg0-1)
               ptconc(idx2) = ptconc(idx2) + 
     &                                   delhg2 * ptconc(idx0) / sumhg0
               ptconc(idx2) = MAX(BNDLPT,ptconc(idx2))
               ptconc(idx0) = ptconc(idx0) - 
     &                                   delhg2 * ptconc(idx2) / sumhg2
               ptconc(idx0) = MAX(BNDLPT,ptconc(idx0))
            enddo
c
c  --- HG2 beign reduced to HG0 ----
c
         else
            ihg0 = iptcls(idxipt(ITRHG0)) - 1
            do i=iptcls(idxipt(ITRHG2)),nptcls(idxipt(ITRHG2))
               idx2 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                      DBLE(numlay) * DBLE(i-1)
               ihg0 = ihg0 + 1
               idx0 = idxcel + DBLE(numcol) * DBLE(numrow) * 
     &                                   DBLE(numlay) * DBLE(ihg0-1)
               ptconc(idx2) = ptconc(idx2) + 
     &                                  delhg2 * ptconc(idx2) / sumhg2
               ptconc(idx2) = MAX(BNDLPT,ptconc(idx2))
               ptconc(idx0) = ptconc(idx0) - 
     &                                  delhg2 * ptconc(idx0) / sumhg0
               ptconc(idx0) = MAX(BNDLPT,ptconc(idx0))
            enddo
         endif
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
