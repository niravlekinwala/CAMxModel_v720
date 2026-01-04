      subroutine soachem(iout,tempk,press,dthr,
     &                   avgox,cprec,cCG,cNV,yfac,
     &                   nfam,avgsen,sen_prec,sen_CG,sen_NV,lddm)
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine solves the oxidation reactions for SOA precursors
c     assuming first-order decay and constant oxidant concentrations
c
c
c     Rate constant:
c        k = A exp(-E/T)                     ; A [cm3/molecule-sec] @ E [K]
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Input arguments:
c        iout     standard output file unit
c        tempk    cell temperature [K]
c        press    cell pressure [mb]
c        dthr     time duration to be integrated [hr]
c        avgox    average oxidant concentrations [ppm]
c                 1 - O; 2 - OH; 3 - O3; 4 - NO3; 5 - NO; 6 - HO2
c        cprec    precursor concentrations [ppm]
c                 1 - BNZA; 2 - TOLA; 3 - XYLA; 4 - IVOA
c                 5 - ISP ; 6 - TRP ; 7 - SQT
c        cCG      CG species concentrations [ppm]
c                 1 - CG1; 2 - CG2; 3 - CG3; 4 - CG4
c        cNV      NV products concentration [ppm]
c                 1 - CGA (anthro); 2 - CGB (bio)
c        yfac     weighted yield factors
c        nfam     number of sensitivity families
c        avgsen   average sensitivities of oxidants [ppm]
c        sen_prec sensitivities of precursors [ppm]
c        sen_CG   sensitivities of CGs [ppm]
c        sen_NV   sensitivities of NVs [ppm]
c        lddm     logical flag for DDM sensitivities
c
c     Output arguments:
c        cprec    precursor concentrations [ppm]
c        cCG      CG species concentrations [ppm]
c        cNV      NV products concentration [ppm]
c        yfac     weighted yield factors
c        sen_prec sensitivities of precursors [ppm]
c        sen_CG   sensitivities of CGs [ppm]
c        sen_NV   sensitivities of NVs [ppm]
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     03/22/06  --bkoo--  Original development
c     06/20/06  --bkoo--  2-product model for isoprene (Henze & Seinfeld, 2006)
c     08/23/13  --bkoo--  Added code for DDM-PM
c     03/18/14  --bkoo--  Revised for benzene SOA
c     08/25/16  --bkoo--  Updated for new SOAP
c     04/22/19  --gy----  Updated hi/low NOx rate constants
c
c-----------------------------------------------------------------------
c
      implicit none
      include 'soap.inc'

      integer :: iout
      real    :: tempk,press,dthr
      real    :: avgox(NOXID),cprec(NPREC_S),cCG(NSOAP),cNV(NNVOL)
      real    :: yfac(3,NPREC_S)

      integer, parameter :: nrxn = 23
      integer, parameter :: idO   = 1, idOH  = 2, idO3  = 3,
     &                      idNO3 = 4, idNO  = 5, idHO2 = 6
      integer, parameter :: idBNZA =  1, idTOLA =  2,
     &                      idXYLA =  3, idIVOA =  4,
     &                      idISP  =  5, idTRP  =  6,
     &                      idSQT  =  7

      real, parameter    :: eps = 1.e-12

      real*8, parameter  :: rkpar(2,nrxn) = reshape( ! rate constant parameters
c                       A         E
     &           (/ 2.30d-12,  190.0d0,    ! rxn  1  BNZA + OH (CB6)
     &              1.80d-12, -340.0d0,    ! rxn  2  TOLA + OH (CB6)
     &              1.85d-11,    0.0d0,    ! rxn  3  XYLA + OH (CB6)
     &              2.70d-12, -360.0d0,    ! rxn  4  ARO-O2 + NO (MCM v3.2 toluene)
     &              2.39d-13,-1300.0d0,    ! rxn  5  ARO-O2 + HO2 (MCM v3.2 toluene)
c
     &              2.70d-11, -390.0d0,    ! rxn  6  ISP + OH (CB6)
     &              1.03d-14, 1995.0d0,    ! rxn  7  ISP + O3 (CB6)
     &              3.03d-12,  448.0d0,    ! rxn  8  ISP + NO3 (CB6)
     &              2.39d-12, -365.0d0,    ! rxn  9  ISP-O2 + NO (CB6)
     &              7.43d-13, -700.0d0,    ! rxn 10  ISP-O2 + HO2 (CB6)
c
     &              1.50d-11, -449.0d0,    ! rxn 11  TRP + OH (CB6)
     &              1.20d-15,  821.0d0,    ! rxn 12  TRP + O3 (CB6)
     &              3.70d-12, -175.0d0,    ! rxn 13  TRP + NO3 (CB6)
     &              2.70d-12, -360.0d0,    ! rxn 14  TRP-O2 + NO (MCM v3.2 a-pinene)
     &              2.66d-13,-1300.0d0,    ! rxn 15  TRP-O2 + HO2 (MCM v3.2 a-pinene)
c
     &              1.97d-10,    0.0d0,    ! rxn 16  SQT + OH (MCM v3.2 b-caryophyllene)
     &              1.20d-14,    0.0d0,    ! rxn 17  SQT + O3 (MCM v3.2 b-caryophyllene)
     &              1.90d-11,    0.0d0,    ! rxn 18  SQT + NO3 (MCM v3.2 b-caryophyllene)
     &              2.70d-12, -360.0d0,    ! rxn 19  SQT-O2 + NO  (assumed same as TRP-O2 + NO)
     &              2.66d-13,-1300.0d0,    ! rxn 20  SQT-O2 + HO2 (assumed same as TRP-O2 + HO2)
c
     &              1.34d-11,    0.0d0,    ! rxn 21  IVOC + OH (Hodzic et al., 2016, ACP)
     &              2.70d-12, -360.0d0,    ! rxn 22  IVOC-O2 + NO  (assumed same as ARO-O2 + NO)
     &              2.39d-13,-1300.0d0 /), ! rxn 23  IVOC-O2 + HO2 (assumed same as ARO-O2 + HO2)
     &           (/ 2, nrxn /) )
c
c======================== DDM Begin =======================
c
      integer :: nfam,ifam
      real    :: avgsen(nfam,NOXID)
      real    :: sen_prec(nfam,NPREC)
      real    :: sen_CG(nfam,NCG)
      real    :: sen_NV(nfam,NNV)
      real    :: srh(nfam)
      real    :: ssumk,srk1(NOXID),sdprec
      logical :: lddm
c
c======================== DDM End =======================
c
c     [1/ppm-hr] = [cm3/molecule-sec]*3600*(6.022e23*1.e-12/8.314)*(P[Pa]/T[K])
      real*8, parameter  :: cf0 = 2.608d14 ! 3600*(6.022e23*1.e-12/8.314)
      real*8  :: cfac

      real    :: rktmp(nrxn), rk1(NOXID), rh, rl, sumk, dprec, cfrac
      real    :: dfac
      integer :: i, j, k

      logical, save :: lfirsttime = .true.
c
c-----Check consistency
c
!      if (lfirsttime) then
!        lfirsttime = .false.
!      endif
c
c-----Conversion factor for rate constants
c
      cfac = cf0 * DBLE(press*100./tempk) ! [1/ppm-hr] = cfac * [cm3/molecule-sec]
c
c-----Rate constants in [1/ppm-hr]
c
      do i = 1, nrxn
        rktmp(i) = SNGL( rkpar(1,i)*dexp(-rkpar(2,i)/DBLE(tempk))*cfac )
      enddo
c
c-----Re-set NV product concentrations to zero
c
      cNV = 0.0
c
c-----NOx-dependent pathways for aromatics
c
      rh = rktmp(4) * avgox(idNO ) ! ARO-O2 + NO
      rl = rktmp(5) * avgox(idHO2) ! ARO-O2 + HO2
      rh = rh / (rh + rl)
      rl = 1.0 - rh

      yfac(1,idBNZA) = y_h(1,idBNZA)*rh + y_l(1,idBNZA)*rl
      yfac(2,idBNZA) = y_h(2,idBNZA)*rh + y_l(2,idBNZA)*rl
      yfac(3,idBNZA) = y_h(3,idBNZA)*rh + y_l(3,idBNZA)*rl

      yfac(1,idTOLA) = y_h(1,idTOLA)*rh + y_l(1,idTOLA)*rl
      yfac(2,idTOLA) = y_h(2,idTOLA)*rh + y_l(2,idTOLA)*rl
      yfac(3,idTOLA) = y_h(3,idTOLA)*rh + y_l(3,idTOLA)*rl

      yfac(1,idXYLA) = y_h(1,idXYLA)*rh + y_l(1,idXYLA)*rl
      yfac(2,idXYLA) = y_h(2,idXYLA)*rh + y_l(2,idXYLA)*rl
      yfac(3,idXYLA) = y_h(3,idXYLA)*rh + y_l(3,idXYLA)*rl
c
c======================== DDM Begin =======================
c
      if ( lddm ) then
        do ifam = 1, nfam
          sen_NV(ifam,:) = 0.0 ! <- cNV reset

          srh(ifam) = ( rl * rktmp(4) * avgsen(ifam,idNO)
     &                - rh * rktmp(5) * avgsen(ifam,idHO2) )
     &                / ( rktmp(4) * avgox(idNO)
     &                  + rktmp(5) * avgox(idHO2) )
c         (rl)' = -(rh)'
        enddo
      endif
c
c======================== DDM End =======================
c
      rk1(idOH ) = rktmp(1) * avgox(idOH ) ! BNZA + OH
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idBNZA) * cfrac
      if (dprec.gt.eps) then
        cprec(idBNZA) = cprec(idBNZA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idBNZA)
        cCG(2) = cCG(2) + dprec * yfac(2,idBNZA)
        cNV(1) = cNV(1) + dprec * yfac(3,idBNZA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rktmp(1) * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idBNZA) * cfrac
     &             + cprec(idBNZA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idBNZA) = sen_prec(ifam,idBNZA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idBNZA)
     &                                      + dprec *
     &            ( y_h(1,idBNZA)*srh(ifam) - y_l(1,idBNZA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idBNZA)
     &                                      + dprec *
     &            ( y_h(2,idBNZA)*srh(ifam) - y_l(2,idBNZA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idBNZA)
     &                                      + dprec *
     &            ( y_h(3,idBNZA)*srh(ifam) - y_l(3,idBNZA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif

      rk1(idOH ) = rktmp(2) * avgox(idOH ) ! TOLA + OH
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idTOLA) * cfrac
      if (dprec.gt.eps) then
        cprec(idTOLA) = cprec(idTOLA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idTOLA)
        cCG(2) = cCG(2) + dprec * yfac(2,idTOLA)
        cNV(1) = cNV(1) + dprec * yfac(3,idTOLA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rktmp(2) * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idTOLA) * cfrac
     &             + cprec(idTOLA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idTOLA) = sen_prec(ifam,idTOLA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idTOLA)
     &                                      + dprec *
     &            ( y_h(1,idTOLA)*srh(ifam) - y_l(1,idTOLA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idTOLA)
     &                                      + dprec *
     &            ( y_h(2,idTOLA)*srh(ifam) - y_l(2,idTOLA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idTOLA)
     &                                      + dprec *
     &            ( y_h(3,idTOLA)*srh(ifam) - y_l(3,idTOLA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif

      rk1(idOH ) = rktmp(3) * avgox(idOH ) ! XYLA + OH
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idXYLA) * cfrac
      if (dprec.gt.eps) then
        cprec(idXYLA) = cprec(idXYLA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idXYLA)
        cCG(2) = cCG(2) + dprec * yfac(2,idXYLA)
        cNV(1) = cNV(1) + dprec * yfac(3,idXYLA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rktmp(3) * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idXYLA) * cfrac
     &             + cprec(idXYLA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idXYLA) = sen_prec(ifam,idXYLA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idXYLA)
     &                                      + dprec *
     &            ( y_h(1,idXYLA)*srh(ifam) - y_l(1,idXYLA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idXYLA)
     &                                      + dprec *
     &            ( y_h(2,idXYLA)*srh(ifam) - y_l(2,idXYLA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idXYLA)
     &                                      + dprec *
     &            ( y_h(3,idXYLA)*srh(ifam) - y_l(3,idXYLA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif

c
c-----NOx-dependent pathways for isoprene
c
      rh = rktmp( 9) * avgox(idNO ) ! ISP-O2 + NO
      rl = rktmp(10) * avgox(idHO2) ! ISP-O2 + HO2
      rh = rh / (rh + rl)
      rl = 1.0 - rh

      yfac(1,idISP) = y_h(1,idISP)*rh + y_l(1,idISP)*rl
      yfac(2,idISP) = y_h(2,idISP)*rh + y_l(2,idISP)*rl
      yfac(3,idISP) = y_h(3,idISP)*rh + y_l(3,idISP)*rl

      rk1(idOH ) = rktmp(6) * avgox(idOH ) ! ISP + OH
      rk1(idO3 ) = rktmp(7) * avgox(idO3 ) ! ISP + O3
      rk1(idNO3) = rktmp(8) * avgox(idNO3) ! ISP + NO3
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)
      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idISP) * cfrac
      if (dprec.gt.eps) then
        cprec(idISP) = cprec(idISP) - dprec
        cCG(3) = cCG(3) + dprec * yfac(1,idISP)
        cCG(4) = cCG(4) + dprec * yfac(2,idISP)
        cNV(2) = cNV(2) + dprec * yfac(3,idISP)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srh(ifam) = ( rl * rktmp( 9) * avgsen(ifam,idNO)
     &                  - rh * rktmp(10) * avgsen(ifam,idHO2) )
     &                  / ( rktmp( 9) * avgox(idNO)
     &                    + rktmp(10) * avgox(idHO2) )
            ssumk = rktmp(6) * avgsen(ifam,idOH)
     &            + rktmp(7) * avgsen(ifam,idO3)
     &            + rktmp(8) * avgsen(ifam,idNO3)
            sdprec = sen_prec(ifam,idISP) * cfrac
     &             + cprec(idISP) * dthr * ssumk ! because cprec has been updated
            sen_prec(ifam,idISP) = sen_prec(ifam,idISP) - sdprec
            sen_CG(ifam,3) = sen_CG(ifam,3) + sdprec * yfac(1,idISP)
     &                                      + dprec *
     &             ( y_h(1,idISP)*srh(ifam) - y_l(1,idISP)*srh(ifam) )
            sen_CG(ifam,4) = sen_CG(ifam,4) + sdprec * yfac(2,idISP)
     &                                      + dprec *
     &             ( y_h(2,idISP)*srh(ifam) - y_l(2,idISP)*srh(ifam) )
            sen_NV(ifam,2) = sen_NV(ifam,2) + sdprec * yfac(3,idISP)
     &                                      + dprec *
     &             ( y_h(3,idISP)*srh(ifam) - y_l(3,idISP)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----NOx-dependent pathways for monoterpenes
c
      rh = rktmp(14) * avgox(idNO ) ! TRP-O2 + NO
      rl = rktmp(15) * avgox(idHO2) ! TRP-O2 + HO2
      rh = rh / (rh + rl)
      rl = 1.0 - rh

      yfac(1,idTRP) = y_h(1,idTRP)*rh + y_l(1,idTRP)*rl
      yfac(2,idTRP) = y_h(2,idTRP)*rh + y_l(2,idTRP)*rl
      yfac(3,idTRP) = y_h(3,idTRP)*rh + y_l(3,idTRP)*rl

      rk1(idOH ) = rktmp(11) * avgox(idOH ) ! TRP + OH
      rk1(idO3 ) = rktmp(12) * avgox(idO3 ) ! TRP + O3
      rk1(idNO3) = rktmp(13) * avgox(idNO3) ! TRP + NO3
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)
      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idTRP) * cfrac
      if (dprec.gt.eps) then
        cprec(idTRP) = cprec(idTRP) - dprec
        cCG(3) = cCG(3) + dprec * yfac(1,idTRP)
        cCG(4) = cCG(4) + dprec * yfac(2,idTRP)
        cNV(2) = cNV(2) + dprec * yfac(3,idTRP)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srh(ifam) = ( rl * rktmp(14) * avgsen(ifam,idNO)
     &                  - rh * rktmp(15) * avgsen(ifam,idHO2) )
     &                  / ( rktmp(14) * avgox(idNO)
     &                    + rktmp(15) * avgox(idHO2) )
            ssumk = rktmp(11) * avgsen(ifam,idOH)
     &            + rktmp(12) * avgsen(ifam,idO3)
     &            + rktmp(13) * avgsen(ifam,idNO3)
            sdprec = sen_prec(ifam,idTRP) * cfrac
     &             + cprec(idTRP) * dthr * ssumk ! because cprec has been updated
            sen_prec(ifam,idTRP) = sen_prec(ifam,idTRP) - sdprec
            sen_CG(ifam,3) = sen_CG(ifam,3) + sdprec * yfac(1,idTRP)
     &                                      + dprec *
     &             ( y_h(1,idTRP)*srh(ifam) - y_l(1,idTRP)*srh(ifam) )
            sen_CG(ifam,4) = sen_CG(ifam,4) + sdprec * yfac(2,idTRP)
     &                                      + dprec *
     &             ( y_h(2,idTRP)*srh(ifam) - y_l(2,idTRP)*srh(ifam) )
            sen_NV(ifam,2) = sen_NV(ifam,2) + sdprec * yfac(3,idTRP)
     &                                      + dprec *
     &             ( y_h(3,idTRP)*srh(ifam) - y_l(3,idTRP)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----NOx-dependent pathways for sesquiterpenes
c
      rh = rktmp(19) * avgox(idNO ) ! SQT-O2 + NO
      rl = rktmp(20) * avgox(idHO2) ! SQT-O2 + HO2
      rh = rh / (rh + rl)
      rl = 1.0 - rh

      yfac(1,idSQT) = y_h(1,idSQT)*rh + y_l(1,idSQT)*rl
      yfac(2,idSQT) = y_h(2,idSQT)*rh + y_l(2,idSQT)*rl
      yfac(3,idSQT) = y_h(3,idSQT)*rh + y_l(3,idSQT)*rl

      rk1(idOH ) = rktmp(16) * avgox(idOH ) ! SQT + OH
      rk1(idO3 ) = rktmp(17) * avgox(idO3 ) ! SQT + O3
      rk1(idNO3) = rktmp(18) * avgox(idNO3) ! SQT + NO3
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)
      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idSQT) * cfrac
      if (dprec.gt.eps) then
        cprec(idSQT) = cprec(idSQT) - dprec
        cCG(3) = cCG(3) + dprec * yfac(1,idSQT)
        cCG(4) = cCG(4) + dprec * yfac(2,idSQT)
        cNV(2) = cNV(2) + dprec * yfac(3,idSQT)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srh(ifam) = ( rl * rktmp(19) * avgsen(ifam,idNO)
     &                  - rh * rktmp(20) * avgsen(ifam,idHO2) )
     &                  / ( rktmp(19) * avgox(idNO)
     &                    + rktmp(20) * avgox(idHO2) )
            ssumk = rktmp(16) * avgsen(ifam,idOH)
     &            + rktmp(17) * avgsen(ifam,idO3)
     &            + rktmp(18) * avgsen(ifam,idNO3)
            sdprec = sen_prec(ifam,idSQT) * cfrac
     &             + cprec(idSQT) * dthr * ssumk ! because cprec has been updated
            sen_prec(ifam,idSQT) = sen_prec(ifam,idSQT) - sdprec
            sen_CG(ifam,3) = sen_CG(ifam,3) + sdprec * yfac(1,idSQT)
     &                                      + dprec *
     &             ( y_h(1,idSQT)*srh(ifam) - y_l(1,idSQT)*srh(ifam) )
            sen_CG(ifam,4) = sen_CG(ifam,4) + sdprec * yfac(2,idSQT)
     &                                      + dprec *
     &             ( y_h(2,idSQT)*srh(ifam) - y_l(2,idSQT)*srh(ifam) )
            sen_NV(ifam,2) = sen_NV(ifam,2) + sdprec * yfac(3,idSQT)
     &                                      + dprec *
     &             ( y_h(3,idSQT)*srh(ifam) - y_l(3,idSQT)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----NOx-dependent pathways for IVOC (IVOA only)
c
      rh = rktmp(22) * avgox(idNO ) ! IVOC-O2 + NO
      rl = rktmp(23) * avgox(idHO2) ! IVOC-O2 + HO2
      rh = rh / (rh + rl)
      rl = 1.0 - rh

      yfac(1,idIVOA) = y_h(1,idIVOA)*rh + y_l(1,idIVOA)*rl
      yfac(2,idIVOA) = y_h(2,idIVOA)*rh + y_l(2,idIVOA)*rl
      yfac(3,idIVOA) = y_h(3,idIVOA)*rh + y_l(3,idIVOA)*rl
c
c======================== DDM Begin =======================
c
      if ( lddm ) then
        do ifam = 1, nfam
          srh(ifam) = ( rl * rktmp(22) * avgsen(ifam,idNO)
     &                - rh * rktmp(23) * avgsen(ifam,idHO2) )
     &                / ( rktmp(22) * avgox(idNO)
     &                  + rktmp(23) * avgox(idHO2) )
        enddo
      endif
c
c======================== DDM End =======================
c
      rk1(idOH ) = rktmp(21) * avgox(idOH ) ! IVOC + OH
      cfrac = 1. - exp(-rk1(idOH )*dthr)

      dprec = cprec(idIVOA) * cfrac
      if (dprec.gt.eps) then
        cprec(idIVOA) = cprec(idIVOA) - dprec
        cCG(1) = cCG(1) + dprec * yfac(1,idIVOA)
        cCG(2) = cCG(2) + dprec * yfac(2,idIVOA)
        cNV(1) = cNV(1) + dprec * yfac(3,idIVOA)
c
c======================== DDM Begin =======================
c
        if ( lddm ) then
          do ifam = 1, nfam
            srk1(idOH) = rktmp(21) * avgsen(ifam,idOH)
            sdprec = sen_prec(ifam,idIVOA) * cfrac
     &             + cprec(idIVOA) * dthr * srk1(idOH) ! because cprec has been updated
            sen_prec(ifam,idIVOA) = sen_prec(ifam,idIVOA) - sdprec
            sen_CG(ifam,1) = sen_CG(ifam,1) + sdprec * yfac(1,idIVOA)
     &                                      + dprec *
     &            ( y_h(1,idIVOA)*srh(ifam) - y_l(1,idIVOA)*srh(ifam) )
            sen_CG(ifam,2) = sen_CG(ifam,2) + sdprec * yfac(2,idIVOA)
     &                                      + dprec *
     &            ( y_h(2,idIVOA)*srh(ifam) - y_l(2,idIVOA)*srh(ifam) )
            sen_NV(ifam,1) = sen_NV(ifam,1) + sdprec * yfac(3,idIVOA)
     &                                      + dprec *
     &            ( y_h(3,idIVOA)*srh(ifam) - y_l(3,idIVOA)*srh(ifam) )
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----Save the weighting factors for NOx-dependent yields for PSAT
c       We assume same rates for ARO-O2 and IVOC-O2
c       If not, PSAT needs to be re-worked
c
      yfac(1,1) = rh
      yfac(2,1) = rl

      return
      end

