      subroutine vbschem(iout,tempk,press,dthr,
     &                   avgox,cprec,cCG,cCGH)
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
c                 1 - BNZA; 2 - TOLA; 3 - XYLA; 4 - ISP; 5 - TRP; 6 - SQT;
c                 7 - IVOG; 8 - IVOD; 9 - IVOA;10 - IVOB
c                 
c        cCG      volatile VBS species concentrations [ppm]
c                 1 - VAS1; 5 - VBS1;  9 - VAP1; 13 - VCP1; 17 - VFP1;
c                 2 - VAS2; 6 - VBS2; 10 - VAP2; 14 - VCP2; 18 - VFP2;
c                 3 - VAS3; 7 - VBS3; 11 - VAP3; 15 - VCP3; 19 - VFP3;
c                 4 - VAS4; 8 - VBS4; 12 - VAP4; 16 - VCP4; 20 - VFP4
c        cCGH     non-volatile VBS species concentrations [ppm]
c                 1 - VAS0; 2 - VBS0;  3 - VAP0;  4 - VCP0;  5 - VFP0
c
c     Output arguments:
c        cprec    precursor concentrations [ppm]
c        cCG      volatile VBS species concentrations [ppm]
c        cCGH     non-volatile VBS species concentrations [ppm]
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     12/07/14  --bkoo--  Modified for VBS from SOACHEM
c     01/08/16  --bkoo--  Updated for Revised VBS
c     08/25/16  --bkoo--  Switched from VBS.INC to SOAP.INC
c     11/28/16  --bkoo--  Revised VBS precursors & reaction rates
c
c-----------------------------------------------------------------------
c
      implicit none
      include 'soap.inc'

      integer :: iout
      real    :: tempk,press,dthr
      real    :: avgox(NOXID),cprec(NPREC_V),cCG(NVBS*NSET),cCGH(NSET)

      integer, parameter :: nrxn = 27

      integer, parameter :: idO   = 1, idOH  = 2, idO3  = 3,
     &                      idNO3 = 4, idNO  = 5, idHO2 = 6

      integer, parameter :: idBNZA =  1, idTOLA =  2,
     &                      idXYLA =  3, idISP  =  4,
     &                      idTRP  =  5, idSQT  =  6,
     &                      idIVOG =  7, idIVOD =  8,
     &                      idIVOA =  9, idIVOB = 10

      real, parameter    :: eps = 1.e-12

      real*8, parameter  :: rkpar(2,nrxn) = reshape( ! rate constant parameters
c                       A         E
     &           (/ 2.30d-12,  190.0d0,    ! rxn  1  BNZA + OH
     &              1.80d-12, -340.0d0,    ! rxn  2  TOLA + OH
     &              1.85d-11,    0.0d0,    ! rxn  3  XYLA + OH
     &              2.70d-12, -360.0d0,    ! rxn  4  ARO-O2 + NO
     &              1.90d-13,-1300.0d0,    ! rxn  5  ARO-O2 + HO2
c
     &              2.70d-11, -390.0d0,    ! rxn  6  ISP + OH
     &              1.03d-14, 1995.0d0,    ! rxn  7  ISP + O3
     &              3.03d-12,  448.0d0,    ! rxn  8  ISP + NO3
     &              2.39d-12, -365.0d0,    ! rxn  9  ISP-O2 + NO
     &              7.43d-13, -700.0d0,    ! rxn 10  ISP-O2 + HO2
c
     &              1.50d-11, -449.0d0,    ! rxn 11  TRP + OH
     &              1.20d-15,  821.0d0,    ! rxn 12  TRP + O3
     &              3.70d-12, -175.0d0,    ! rxn 13  TRP + NO3
     &              8.80d-13, -180.2d0,    ! rxn 14  TRP-O2 + NO
     &              4.10d-13, -790.0d0,    ! rxn 15  TRP-O2 + HO2
c
     &              1.97d-10,    0.0d0,    ! rxn 16  SQT + OH
     &              1.16d-14,    0.0d0,    ! rxn 17  SQT + O3
     &              1.90d-11,    0.0d0,    ! rxn 18  SQT + NO3
c
     &              4.00d-11,    0.0d0,    ! rxn 19  IVOG + OH
     &              4.00d-11,    0.0d0,    ! rxn 20  IVOD + OH
     &              4.00d-11,    0.0d0,    ! rxn 21  IVOA + OH - temporary
     &              4.00d-11,    0.0d0,    ! rxn 22  IVOB + OH
c
     &              2.00d-11,    0.0d0,    ! rxn 23  VAS? + OH
     &                 0.0d0,    0.0d0,    ! rxn 24  VBS? + OH
     &              4.00d-11,    0.0d0,    ! rxn 25  VAP? + OH
     &              4.00d-11,    0.0d0,    ! rxn 26  VCP? + OH
     &              4.00d-11,    0.0d0 /), ! rxn 27  VFP? + OH
     &           (/ 2, nrxn /) )

      real, parameter    :: yh(0:nvbs,NPREC_V) = reshape( ! high-NOx yield [ppm/ppm]
c                    bin0   bin1   bin2   bin3   bin4
     &           (/ 0.0  , 0.001, 0.079, 0.148, 0.222,    ! BNZA
     &              0.0  , 0.006, 0.145, 0.281, 0.432,    ! TOLA
     &              0.0  , 0.001, 0.127, 0.201, 0.301,    ! XYLA
     &              0.0  , 0.0  , 0.009, 0.006, 0.0  ,    ! ISP
     &              0.0  , 0.010, 0.101, 0.173, 0.451,    ! TRP
     &              0.0  , 0.092, 0.188, 0.968, 0.679,    ! SQT  - NOx-independent
     &              0.022, 0.098, 0.373, 0.699, 0.0  ,    ! IVOG - NOx-independent
     &              0.081, 0.135, 0.800, 0.604, 0.0  ,    ! IVOD - NOx-independent
     &              0.081, 0.135, 0.800, 0.604, 0.0  ,    ! IVOA - NOx-independent - temporary
     &              0.081, 0.135, 0.800, 0.604, 0.0   /), ! IVOB - NOx-independent
     &           (/ nvbs+1, NPREC_V /) )

      real, parameter    :: yl(0:nvbs,NPREC_V) = reshape( ! low-NOx yield [ppm/ppm]
c                    bin0   bin1   bin2   bin3   bin4
     &           (/ 0.0  , 0.035, 0.108, 0.185, 0.268,    ! BNZA
     &              0.0  , 0.006, 0.145, 0.437, 0.281,    ! TOLA
     &              0.0  , 0.048, 0.195, 0.252, 0.364,    ! XYLA
     &              0.0  , 0.004, 0.013, 0.006, 0.0  ,    ! ISP
     &              0.0  , 0.087, 0.077, 0.309, 0.540,    ! TRP
     &              0.0  , 0.092, 0.188, 0.968, 0.679,    ! SQT  - NOx-independent
     &              0.022, 0.098, 0.373, 0.699, 0.0  ,    ! IVOG - NOx-independent
     &              0.081, 0.135, 0.800, 0.604, 0.0  ,    ! IVOD - NOx-independent
     &              0.081, 0.135, 0.800, 0.604, 0.0  ,    ! IVOA - NOx-independent - temporary
     &              0.081, 0.135, 0.800, 0.604, 0.0   /), ! IVOB - NOx-independent
     &           (/ nvbs+1, NPREC_V /) )

      real, parameter    :: yy2(0:nvbs) =               ! TRP + NO3 yield [ppm/ppm]
c                    bin0   bin1   bin2   bin3   bin4
     &           (/ 0.314, 0.029, 0.0  , 0.862, 0.0   /)

      real, parameter    :: ya(2,nvbs) = reshape( ! yield for VAP aging
c                     VAS    VAP
     &           (/ 0.142, 0.864,    ! bin0
     &              0.129, 0.877,    ! bin1
     &              0.116, 0.889,    ! bin2
     &              0.137, 0.869 /), ! bin3
     &           (/ 2, nvbs /) )

      real, parameter    :: yf(2,nvbs) = reshape( ! yield for VFP aging
c                     VBS    VFP
     &           (/ 0.464, 0.538,    ! bin0
     &              0.313, 0.689,    ! bin1
     &              0.220, 0.783,    ! bin2
     &              0.156, 0.846 /), ! bin3
     &           (/ 2, nvbs /) )

c     [1/ppm-hr] = [cm3/molecule-sec]*3600*(6.022e23*1.e-12/8.314)*(P[Pa]/T[K])
      real*8, parameter  :: cf0 = 2.608d14 ! 3600*(6.022e23*1.e-12/8.314)
      real*8  :: cfac

      real    :: cVAS(0:nvbs),cVBS(0:nvbs),
     &           cVAP(0:nvbs),cVCP(0:nvbs),cVFP(0:nvbs)
      real    :: rktmp(nrxn), rk1(NOXID), rh, rl, rr, rr2, sumk, dprec
      real    :: fage, cfrac
      integer :: i, j, k, nchk, iprec

      logical, save :: lfirsttime = .true.
c
c-----Check consistency
c
!      if (lfirsttime) then
!        lfirsttime = .false.
!      endif
c
c-----Load VBS species
c
      cVAS(0) = 0.0
      cVBS(0) = 0.0
      cVAP(0) = 0.0
      cVCP(0) = 0.0
      cVFP(0) = 0.0
      do j = 1, nvbs
        cVAS(j) = cCG(j)
        cVBS(j) = cCG(j+nvbs)
        cVAP(j) = cCG(j+2*nvbs)
        cVCP(j) = cCG(j+3*nvbs)
        cVFP(j) = cCG(j+4*nvbs)
      enddo
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
c-----Photochemical againg
c
      ! VAS?
      rk1(idOH ) = rktmp(23) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      do j = 1, nvbs
        dprec = cVAS(j) * cfrac
        if (dprec.gt.eps) then
          cVAS(j-1) = cVAS(j-1) + dprec
          cVAS(j)   = cVAS(j)   - dprec
        endif
      enddo
      ! VBS?
      rk1(idOH ) = rktmp(24) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      do j = 1, nvbs
        dprec = cVBS(j) * cfrac
        if (dprec.gt.eps) then
          cVBS(j-1) = cVBS(j-1) + dprec
          cVBS(j)   = cVBS(j)   - dprec
        endif
      enddo
      ! VAP?
      rk1(idOH ) = rktmp(25) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      do j = 1, nvbs
        dprec = cVAP(j) * cfrac
        if (dprec.gt.eps) then
          cVAS(j-1) = cVAS(j-1) + dprec * ya(1,j)
          cVAP(j-1) = cVAP(j-1) + dprec * ya(2,j)
          cVAP(j)   = cVAP(j)   - dprec
        endif
      enddo
      ! VCP?
      rk1(idOH ) = rktmp(26) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      do j = 1, nvbs
        dprec = cVCP(j) * cfrac
        if (dprec.gt.eps) then
          cVBS(j-1) = cVBS(j-1) + dprec * ya(1,j)
          cVCP(j-1) = cVCP(j-1) + dprec * ya(2,j)
          cVCP(j)   = cVCP(j)   - dprec
        endif
      enddo
      ! VFP?
      rk1(idOH ) = rktmp(27) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      do j = 1, nvbs
        dprec = cVFP(j) * cfrac
        if (dprec.gt.eps) then
          cVBS(j-1) = cVBS(j-1) + dprec * yf(1,j)
          cVFP(j-1) = cVFP(j-1) + dprec * yf(2,j)
          cVFP(j)   = cVFP(j)   - dprec
        endif
      enddo
c
c-----NOx-dependent pathways for aromatics
c
      rh = rktmp(4) * avgox(idNO )
      rl = rktmp(5) * avgox(idHO2)
      rh = rh / (rh + rl)
      rl = 1.0 - rh
      ! BNZA
      rk1(idOH ) = rktmp(1) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idBNZA) * cfrac
      if (dprec.gt.eps) then
        cprec(idBNZA) = cprec(idBNZA) - dprec
        do j = 0, nvbs
          cVAS(j) = cVAS(j) + dprec *
     &              ( yh(j,idBNZA)*rh + yl(j,idBNZA)*rl )
        enddo
      endif
      ! TOLA
      rk1(idOH ) = rktmp(2) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idTOLA) * cfrac
      if (dprec.gt.eps) then
        cprec(idTOLA) = cprec(idTOLA) - dprec
        do j = 0, nvbs
          cVAS(j) = cVAS(j) + dprec *
     &              ( yh(j,idTOLA)*rh + yl(j,idTOLA)*rl )
        enddo
      endif
      ! XYLA
      rk1(idOH ) = rktmp(3) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idXYLA) * cfrac
      if (dprec.gt.eps) then
        cprec(idXYLA) = cprec(idXYLA) - dprec
        do j = 0, nvbs
          cVAS(j) = cVAS(j) + dprec *
     &              ( yh(j,idXYLA)*rh + yl(j,idXYLA)*rl )
        enddo
      endif
c
c-----NOx-dependent pathways for isoprene
c
      rh = rktmp( 9) * avgox(idNO )
      rl = rktmp(10) * avgox(idHO2)
      rh = rh / (rh + rl)
      rl = 1.0 - rh

      rk1(idOH ) = rktmp(6) * avgox(idOH )
      rk1(idO3 ) = rktmp(7) * avgox(idO3 )
      rk1(idNO3) = rktmp(8) * avgox(idNO3)
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)
      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idISP ) * cfrac
      if (dprec.gt.eps) then
        cprec(idISP ) = cprec(idISP ) - dprec
        do j = 0, nvbs
          cVBS(j) = cVBS(j) + dprec *
     &              ( yh(j,idISP )*rh + yl(j,idISP )*rl )
        enddo
      endif
c
c-----NOx-dependent pathways for monoterpenes
c
      rh = rktmp(14) * avgox(idNO )
      rl = rktmp(15) * avgox(idHO2)
      rh = rh / (rh + rl)
      rl = 1.0 - rh

      rk1(idOH ) = rktmp(11) * avgox(idOH )
      rk1(idO3 ) = rktmp(12) * avgox(idO3 )
      rk1(idNO3) = rktmp(13) * avgox(idNO3)
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)
      rr2 = rk1(idNO3) / sumk
      rr = 1. - rr2
      cfrac = 1. - exp(-sumk*dthr)
      dprec = cprec(idTRP ) * cfrac
      if (dprec.gt.eps) then
        cprec(idTRP ) = cprec(idTRP ) - dprec
        do j = 0, nvbs
          cVBS(j) = cVBS(j) + dprec * rr *
     &              ( yh(j,idTRP )*rh + yl(j,idTRP )*rl )
          cVBS(j) = cVBS(j) + dprec * rr2 * yy2(j)
        enddo
      endif
c
c-----NOx-independent pathways for sesquiterpenes
c
      rk1(idOH ) = rktmp(16) * avgox(idOH )
      rk1(idO3 ) = rktmp(17) * avgox(idO3 )
      rk1(idNO3) = rktmp(18) * avgox(idNO3)
      sumk = rk1(idOH ) + rk1(idO3 ) + rk1(idNO3)
      dprec = cprec(idSQT ) * ( 1. - exp(-sumk*dthr) )
      if (dprec.gt.eps) then
        cprec(idSQT ) = cprec(idSQT ) - dprec
        do j = 0, nvbs
          cVBS(j) = cVBS(j) + dprec * yh(j,idSQT )
        enddo
      endif
c
c-----NOx-independent pathways for IVOC
c
      ! IVOG
      rk1(idOH ) = rktmp(19) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idIVOG) * cfrac
      if (dprec.gt.eps) then
        cprec(idIVOG) = cprec(idIVOG) - dprec
        do j = 0, nvbs
          cVAS(j) = cVAS(j) + dprec * yh(j,idIVOG)
        enddo
      endif
      ! IVOD
      rk1(idOH ) = rktmp(20) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idIVOD) * cfrac
      if (dprec.gt.eps) then
        cprec(idIVOD) = cprec(idIVOD) - dprec
        do j = 0, nvbs
          cVAS(j) = cVAS(j) + dprec * yh(j,idIVOD)
        enddo
      endif
      ! IVOA
      rk1(idOH ) = rktmp(21) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idIVOA) * cfrac
      if (dprec.gt.eps) then
        cprec(idIVOA) = cprec(idIVOA) - dprec
        do j = 0, nvbs
          cVAS(j) = cVAS(j) + dprec * yh(j,idIVOA)
        enddo
      endif
      ! IVOB
      rk1(idOH ) = rktmp(22) * avgox(idOH )
      cfrac = 1. - exp(-rk1(idOH )*dthr)
      dprec = cprec(idIVOB) * cfrac
      if (dprec.gt.eps) then
        cprec(idIVOB) = cprec(idIVOB) - dprec
        do j = 0, nvbs
          cVBS(j) = cVBS(j) + dprec * yh(j,idIVOB)
        enddo
      endif
c
c-----Put VBS species back
c
      cCGH(1) = cVAS(0)
      cCGH(2) = cVBS(0)
      cCGH(3) = cVAP(0)
      cCGH(4) = cVCP(0)
      cCGH(5) = cVFP(0)
      do j = 1, nvbs
        cCG(j)        = cVAS(j)
        cCG(j+nvbs)   = cVBS(j)
        cCG(j+2*nvbs) = cVAP(j)
        cCG(j+3*nvbs) = cVCP(j)
        cCG(j+4*nvbs) = cVFP(j)
      enddo

      return
      end

