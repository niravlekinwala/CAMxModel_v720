C=====================================================================
C
C   DDM-PM SUBROUTINES FOR RADM
C
C   Subroutines included in this file:
C     DDMRADINI
C     DDMRADM
C     DDMRADUNIT
C     DDMRADEQL
C     DDMRADSUB1
C     DDMRADSUB2
C     DDMDAVIES
C     DDMRADRXN
C
C   Dependent upon:
C     tracer.com -> camx.prm  - for MXFDDM, NDDMSP & LDDM
C     ddmradm.inc
C
C   Log:
C     04/26/2005  bkoo - original development
C     06/25/2008  bkoo - added optional model species (FOA/MHP/OHP/PAA/OPA)
C     09/19/2013  bkoo - revised for CaCO3 which is now a model species
C     10/26/2018  bkoo - updated for revised RADM
C
C=====================================================================

      SUBROUTINE DDMRADINI (CSO2,CNA,CCL,fCA_PRM,fCA_CRS,CONVFAC,
     &                      WTMOL,SDDM,NFAM,NSEN,NWTMOL)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmradm.inc'

      INTEGER   NFAM, NSEN, NWTMOL
      REAL      SDDM(NFAM,NSEN)           ! [ppm]/[p] for gas; [ug/m3]/[p] for PM
      REAL      WTMOL(NWTMOL)             ! Molecular weight array
      REAL      CONVFAC                   ! umol/m3 = ppm * convfac
      REAL      fCA_PRM, fCA_CRS          ! mass fraction of Ca on FPRM/FCRS
      REAL      CNA, CCL                  ! CNA = sodium conc; CCL = chloride conc
      REAL      CSO2                      ! SO2 [ppm]

      INTEGER   J
C
C     Check if SO2 is too small
C
      IF ( CSO2.LT.1.E-6 ) THEN
        IPAR = 0
        RETURN
      ENDIF
C
C     Convert to [mol/mol of air]/[p]
C
      DO J = 1, NFAM
C     LOADDMPM has initialized SDDM rows for non-model species to zero
        SINI(jTNH3 ,J) = ( SDDM(J,jdpNH3 ) +
     &                     SDDM(J,jdpPNH4)/WTMOL(2)/CONVFAC )*1.E-6
        SINI(jTHNO3,J) = ( SDDM(J,jdpHNO3) + SDDM(J,jdpN2O5)*2.0 +
     &                     SDDM(J,jdpPNO3)/WTMOL(3)/CONVFAC )*1.E-6
        SINI(jTFOA ,J) =   SDDM(J,jdpFOA )*1.E-6
        SINI(jTCO2 ,J) = 0.0
        SINI(jTSO2 ,J) =   SDDM(J,jdpSO2 )*1.E-6
        SINI(jTSVI ,J) = ( SDDM(J,jdpSULF) +
     &                     SDDM(J,jdpPSO4)/WTMOL(1)/CONVFAC )*1.E-6
        SINI(jTH2O2,J) =   SDDM(J,jdpH2O2)*1.E-6
        SINI(jTO3  ,J) =   SDDM(J,jdpO3  )*1.E-6
        SINI(jTMHP ,J) = ( SDDM(J,jdpMHP ) + SDDM(J,jdpOHP ) )*1.E-6
        SINI(jTPAA ,J) = ( SDDM(J,jdpPAA ) + SDDM(J,jdpOPA ) )*1.E-6
        SINI(jTGLY ,J) =   SDDM(J,jdpGLY )*1.E-6
        SINI(jTMGLY,J) =   SDDM(J,jdpMGLY)*1.E-6
        SINI(jTGLYD,J) =   SDDM(J,jdpGLYD)*1.E-6
        SINI(jTHO  ,J) =   SDDM(J,jdpHO  )*1.E-6
        IF (CNA.GT.CCL) THEN
          SINI(jTNACL,J) = SDDM(J,jdpPCL )/WTMOL(4)/CONVFAC*1.E-6
        ELSE
          SINI(jTNACL,J) = SDDM(J,jdpNA  )/WTMOL(5)/CONVFAC*1.E-6
        ENDIF
        ! assumed Ca is scaled from FPRM+FCRS
        ! Fe, Mn, Mg & K are assumed background constants
        SINI(jTCAC ,J) = (fCA_PRM*SDDM(J,jdpFPRM)+
     &                    fCA_CRS*SDDM(J,jdpFCRS))/WTMOL(8)/CONVFAC*1.E-6
        SINI(jTSOAC,J) = SDDM(J,jdpSOPB)/WTMOL(12)/CONVFAC*1.E-6

        SAVSISO2(J) = SINI(jTSO2,J)
        SAVSISVI(J) = SINI(jTSVI,J)

      ENDDO

      IPAR = 1

      RETURN
      END
C     END OF SUBROUTINE DDMRADINI


      SUBROUTINE DDMRADM (r_1,r_11,r_2,r_22,
     &                    CONVFAC,WTMOL,SDDM,NFAM,NSEN,NWTMOL,IERR,IOUT)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmradm.inc'

      INTEGER   NFAM, NSEN, NWTMOL
      REAL      SDDM(NFAM,NSEN)           ! [ppm]/[p] for gas; [ug/m3]/[p] for PM
      REAL      WTMOL(NWTMOL)             ! Molecular weight array
      REAL      CONVFAC                   ! umol/m3 = ppm * convfac
      REAL      r_1                       ! con(MHP    ;t)   /con(MHP+OHP;t)
      REAL      r_11                      ! con(MHP+OHP;t+dt)/con(MHP+OHP;t)
      REAL      r_2                       ! con(PAA    ;t)   /con(PAA+OPA;t)
      REAL      r_22                      ! con(PAA+OPA;t+dt)/con(PAA+OPA;t)
      INTEGER   IERR, IOUT

      INTEGER   J
C
C     Initialize error flag
C
      IERR = 1
C
C     Check if SO2 is too small
C
      IF ( IPAR.EQ.0 ) THEN
        DO J = 1, NFAM
          SDDM(J,jdpPSO4) = SDDM(J,jdpPSO4)
     &                    + SDDM(J,jdpSULF)*WTMOL(1)*CONVFAC
          SDDM(J,jdpHNO3) = SDDM(J,jdpHNO3) + SDDM(J,jdpN2O5)*2.0
          SDDM(J,jdpSULF) = 0.0
          SDDM(J,jdpN2O5) = 0.0
        ENDDO
        IERR = IPAR
        RETURN
      ENDIF
C
C     Check error from the previous routines
C
      IF ( IPAR.LT.0 ) THEN
        IERR = IPAR
        RETURN
      ENDIF
C
C     Assign the solution
C
      DO J = 1, NFAM
C     Currently sensitivities of HNO3/NO3 and NH3/NH4 are not properly
C     implemented since RADM uses initial partition factors to partition
C     them and total (gas+aerosol) does not change. Also, they will be
C     re-partitioned by the following ISORROPIA anyway. However, it should
C     be considered here that all N2O5 is added to HNO3/NO3 by RADM.
        SDDM(J,jdpHNO3) = SDDM(J,jdpHNO3) + SDDM(J,jdpN2O5)*2.0
        SDDM(J,jdpN2O5) = 0.0
        SDDM(J,jdpSULF) = 0.0
        SDDM(J,jdpPSO4) = SINI(jTSVI ,J)*WTMOL(1)*CONVFAC*1.E6
        SDDM(J,jdpSO2 ) = SINI(jTSO2 ,J)*1.E6
        SDDM(J,jdpH2O2) = SINI(jTH2O2,J)*1.E6
        SDDM(J,jdpO3  ) = SINI(jTO3  ,J)*1.E6
C     For CB05 and SAPRC99, FOA/MHP/OHP/PAA/OPA are (optional) model species.
C     However, FOA is not changed by RADM, thus neither is the sensitivity.
C     CaCO3 is now also a model species, but not changed by RADM.
        SDDM(J,jdpMHP ) = SINI(jTMHP ,J)*1.E6*r_1 + r_11*
     &                   (SDDM(J,jdpMHP )*(1.-r_1)-SDDM(J,jdpOHP )*r_1)
        SDDM(J,jdpOHP ) = SINI(jTMHP ,J)*1.E6 - SDDM(J,jdpMHP )
        SDDM(J,jdpPAA ) = SINI(jTPAA ,J)*1.E6*r_2 + r_22*
     &                   (SDDM(J,jdpPAA )*(1.-r_2)-SDDM(J,jdpOPA )*r_2)
        SDDM(J,jdpOPA ) = SINI(jTPAA ,J)*1.E6 - SDDM(J,jdpPAA )
        SDDM(J,jdpGLY ) = SINI(jTGLY ,J)*1.E6
        SDDM(J,jdpMGLY) = SINI(jTMGLY,J)*1.E6
        SDDM(J,jdpGLYD) = SINI(jTGLYD,J)*1.E6
        ! RADM does not consume OH radical
        SDDM(J,jdpSOPB) = SINI(jTSOAC,J)*WTMOL(12)*CONVFAC*1.E6
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMRADM


      SUBROUTINE DDMRADUNIT (IEND,PRESATMOVERXL)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmradm.inc'

      INTEGER   IEND          ! 0 - at the start; 1 - at the end
      REAL      PRESATMOVERXL ! [M] = PRESATMOVERXL * [mol/mol of air]
      REAL      CONVFAC
      INTEGER   I, J
C
C     Initial check
C
      IF ( .NOT.LDDM ) RETURN             ! Check if DDM is being used
      IF ( IPAR.LE.0 ) RETURN             ! Check error from DDMRADEQL
C
C     Set conversion factor
C
      CONVFAC = PRESATMOVERXL
      IF (IEND.EQ.1) CONVFAC = 1. / CONVFAC
C
C     Unit conversion: [mol/mol of air]/[p] <-> [M]/[p]
C
      DO J = 1, NDDMSP
        DO I = 1, NINPT
          SINI(I,J) = SINI(I,J) * CONVFAC
        ENDDO
        SAVSISO2(J) = SAVSISO2(J) * CONVFAC
        SAVSISVI(J) = SAVSISVI(J) * CONVFAC
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMRADUNIT


      SUBROUTINE DDMRADEQL (NLIQS,LIQUID,STION,GM1,GM2,
     &                      XLNH3,XLHNO3,XLFOA,XLCO2,XLSO2,
     &                      XLH2O2,XLO3,XLMHP,XLPAA,
     &                      XLGLY,XLMGLY,XLGLYD,XLHO,IOUT)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmradm.inc'

      INTEGER   NLIQS
      REAL      LIQUID(NLIQS) ! Liquid concentration array [M]
      REAL      STION         ! Ionic strength
      REAL      GM1           ! Activity coefficient for ions with absolute charge = 1
      REAL      GM2           ! Activity coefficient for ions with absolute charge = 2
      REAL      XLNH3         ! = NH3H  * XL
      REAL      XLHNO3        ! = HNO3H * XL
      REAL      XLFOA         ! = FOAH  * XL (new)
      REAL      XLCO2         ! = CO2H  * XL
      REAL      XLSO2         ! = SO2H  * XL
      REAL      XLH2O2        ! = H2O2H * XL
      REAL      XLO3          ! = O3H   * XL
      REAL      XLMHP         ! = MHPH  * XL
      REAL      XLPAA         ! = PAAH  * XL
      REAL      XLGLY         ! = GLYH  * XL
      REAL      XLMGLY        ! = MGLYH * XL
      REAL      XLGLYD        ! = GLYDH * XL
      REAL      XLHO          ! = HOH   * XL
      INTEGER   IOUT          ! Unit number of log file

      LOGICAL   LSEQ(NSENS)               ! If true, the eqn is solved
      LOGICAL   LSEN(NSENS)               ! If true, the sensitivity is solved
      REAL      AMAT(NSENS,NSENS)         ! Left-hand side matrix A
                                          ! Right-hand side vector B -> XMAT
      REAL      AA(NSENS,NSENS),BB(NSENS) ! These are passed to DDMEQNSLV solely
      INTEGER   IPVT(NSENS)               ! to allow their size to be adjustable

      REAL      SION, DGM11, DGM12, DGM21, DGM22
      REAL      CX, C1
      INTEGER   NDIM, iEQ, I, J, K
C
C     Initial check
C
      IF ( .NOT.LDDM ) RETURN             ! Check if DDM is being used
      IF ( IPAR.LE.0 ) RETURN             ! Check error from the previous call
C
C     Initialize arrays
C
      Y(jH     ) = LIQUID( 1)
      Y(jNH4   ) = LIQUID( 2)
      Y(jNO3   ) = LIQUID(14)
      Y(jHCO2  ) = LIQUID(23)
      Y(jHCO3  ) = LIQUID(12)
      Y(jHSO3  ) = LIQUID( 9)
      Y(jHSO4  ) = LIQUID( 7)
      Y(jOH    ) = LIQUID( 5)
      Y(jCO3   ) = LIQUID(11)
      Y(jSO3   ) = LIQUID( 8)
      Y(jSO4   ) = LIQUID( 6)
      Y(jNH3aq ) = LIQUID(15)
      Y(jHNO3aq) = LIQUID(32)
      Y(jFOAaq ) = LIQUID(22)
      Y(jCO2aq ) = LIQUID(13)
      Y(jSO2aq ) = LIQUID(10)
      Y(jH2O2aq) = LIQUID(17)
      Y(jO3aq  ) = LIQUID(18)
      Y(jMHPaq ) = LIQUID(24)
      Y(jPAAaq ) = LIQUID(25)
      Y(jGLYaq ) = LIQUID(34)
      Y(jMGLYaq) = LIQUID(35)
      Y(jGLYDaq) = LIQUID(36)
      Y(jHOaq  ) = LIQUID(37)
      DO I = 1, NSENS
        DO J = 1, NSENS
          AMAT(J,I) = 0.0
        ENDDO
      ENDDO
C     Initialize the solution matrix
C     - should be done here as this routine is called every (internal) timestep
      DO J = 1, NDDMSP
        DO I = 1, NSENS
          XMAT(I,J) = 0.0
        ENDDO
      ENDDO
      NDIM = 16 ! Should be consistent with the number of equations
      DO I = 1, NDIM
        LSEQ(I) = .TRUE.
        LSEN(I) = .TRUE.
      ENDDO
      DO I = NDIM + 1, NSENS
        LSEQ(I) = .FALSE.
        LSEN(I) = .FALSE.
      ENDDO
C
C     Select eqns & sensitivities to solve
C
      IF ( Y(jSO2aq ).LE.0.0 ) THEN
        CALL DDMEQNSUB (LSEQ(iMBSO2 ),LSEN(jSO2aq ))
        CALL DDMEQNSUB (LSEQ(iKS1   ),LSEN(jHSO3  ))
        CALL DDMEQNSUB (LSEQ(iKS2   ),LSEN(jSO3   ))
        NDIM = NDIM - 3
        DO J = 1, NDDMSP
          SINI(jTSO2,J) = 0.0 ! TSO2 is depleted
          SINI(jTSVI,J) = SAVSISVI(J) + SAVSISO2(J)
        ENDDO
      ENDIF
      IF ( Y(jNH3aq ).LE.0.0 ) THEN
        CALL DDMEQNSUB (LSEQ(iMBNH3 ),LSEN(jNH3aq ))
        CALL DDMEQNSUB (LSEQ(iKA    ),LSEN(jNH4   ))
        NDIM = NDIM - 2
      ENDIF
      IF ( Y(jHNO3aq).LE.0.0 ) THEN
        CALL DDMEQNSUB (LSEQ(iMBHNO3),LSEN(jHNO3aq))
        CALL DDMEQNSUB (LSEQ(iKN    ),LSEN(jNO3   ))
        NDIM = NDIM - 2
      ENDIF
      IF ( Y(jFOAaq ).LE.0.0 ) THEN
        CALL DDMEQNSUB (LSEQ(iMBFOA ),LSEN(jFOAaq ))
        CALL DDMEQNSUB (LSEQ(iKF    ),LSEN(jHCO2  ))
        NDIM = NDIM - 2
      ENDIF
C====================================================================
C     Fill the matrix A and constant vector B
C====================================================================
C
C     Activity coefficients
C
      SION = SQRT(STION)
      CALL DDMDAVIES (SION,GM1,GM2,DGM11,DGM12,DGM21,DGM22)
      CX = 1.     ! Equlibrium equations are now normalized
C
C     [H2O(aq)] <-> [H+] + [OH-]
C
      iEQ = iKW
CCC   IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jOH)*GM1*GM1
        C1 = 2.*CX/GM1
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM11,DGM12,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM11)      ! Constant term
        AMAT(jH ,iEQ) = AMAT(jH ,iEQ) + CX/Y(jH)
        AMAT(jOH,iEQ) = AMAT(jOH,iEQ) + CX/Y(jOH)
CCC   ENDIF
C
C     [NH3(aq)] <-> [NH4+] + [OH-]
C
      iEQ = iKA
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jNH4)*Y(jOH)*GM1*GM1/Y(jNH3aq)
        C1 = 2.*CX/GM1
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM11,DGM12,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM11)      ! Constant term
        AMAT(jNH4  ,iEQ) = AMAT(jNH4  ,iEQ) + CX/Y(jNH4)
        AMAT(jOH   ,iEQ) = AMAT(jOH   ,iEQ) + CX/Y(jOH)
        AMAT(jNH3aq,iEQ) = AMAT(jNH3aq,iEQ) - CX/Y(jNH3aq)
      ENDIF
C
C     [HNO3(aq)] <-> [H+] + [NO3-]
C
      iEQ = iKN
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jNO3)*GM1*GM1/Y(jHNO3aq)
        C1 = 2.*CX/GM1
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM11,DGM12,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM11)      ! Constant term
        AMAT(jH     ,iEQ) = AMAT(jH     ,iEQ) + CX/Y(jH)
        AMAT(jNO3   ,iEQ) = AMAT(jNO3   ,iEQ) + CX/Y(jNO3)
        AMAT(jHNO3aq,iEQ) = AMAT(jHNO3aq,iEQ) - CX/Y(jHNO3aq)
      ENDIF
C
C     [HCOOH(aq)] <-> [H+] + [HCOO-]
C
      iEQ = iKF
CCC   IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jHCO2)*GM1*GM1/Y(jFOAaq)
        C1 = 2.*CX/GM1
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM11,DGM12,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM11)      ! Constant term
        AMAT(jH    ,iEQ) = AMAT(jH    ,iEQ) + CX/Y(jH)
        AMAT(jHCO2 ,iEQ) = AMAT(jHCO2 ,iEQ) + CX/Y(jHCO2)
        AMAT(jFOAaq,iEQ) = AMAT(jFOAaq,iEQ) - CX/Y(jFOAaq)
CCC   ENDIF
C
C     [CO2(aq)] <-> [H+] + [HCO3-]
C
      iEQ = iKC1
CCC   IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jHCO3)*GM1*GM1/Y(jCO2aq)
        C1 = 2.*CX/GM1
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM11,DGM12,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM11)      ! Constant term
        AMAT(jH    ,iEQ) = AMAT(jH    ,iEQ) + CX/Y(jH)
        AMAT(jHCO3 ,iEQ) = AMAT(jHCO3 ,iEQ) + CX/Y(jHCO3)
        AMAT(jCO2aq,iEQ) = AMAT(jCO2aq,iEQ) - CX/Y(jCO2aq)
CCC   ENDIF
C
C     [HCO3-] <-> [H+] + [CO3=]
C
      iEQ = iKC2
CCC   IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jCO3)*GM2/Y(jHCO3)
        C1 = CX/GM2
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM21,DGM22,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM21)      ! Constant term
        AMAT(jH   ,iEQ) = AMAT(jH   ,iEQ) + CX/Y(jH)
        AMAT(jCO3 ,iEQ) = AMAT(jCO3 ,iEQ) + CX/Y(jCO3)
        AMAT(jHCO3,iEQ) = AMAT(jHCO3,iEQ) - CX/Y(jHCO3)
CCC   ENDIF
C
C     [SO2(aq)] <-> [H+] + [HSO3-]
C
      iEQ = iKS1
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jHSO3)*GM1*GM1/Y(jSO2aq)
        C1 = 2.*CX/GM1
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM11,DGM12,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM11)      ! Constant term
        AMAT(jH    ,iEQ) = AMAT(jH    ,iEQ) + CX/Y(jH)
        AMAT(jHSO3 ,iEQ) = AMAT(jHSO3 ,iEQ) + CX/Y(jHSO3)
        AMAT(jSO2aq,iEQ) = AMAT(jSO2aq,iEQ) - CX/Y(jSO2aq)
      ENDIF
C
C     [HSO3-] <-> [H+] + [SO3=]
C
      iEQ = iKS2
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jSO3)*GM2/Y(jHSO3)
        C1 = CX/GM2
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM21,DGM22,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM21)      ! Constant term
        AMAT(jH   ,iEQ) = AMAT(jH   ,iEQ) + CX/Y(jH)
        AMAT(jSO3 ,iEQ) = AMAT(jSO3 ,iEQ) + CX/Y(jSO3)
        AMAT(jHSO3,iEQ) = AMAT(jHSO3,iEQ) - CX/Y(jHSO3)
      ENDIF
C
C     [HSO4-] <-> [H+] + [SO4=]
C
      iEQ = iKS3
CCC   IF ( LSEQ(iEQ) ) THEN
CX1        CX = Y(jH)*Y(jSO4)*GM2/Y(jHSO4)
        C1 = CX/GM2
        CALL DDMRADSUB1 (NSENS,jOH,NIONS,C1,DGM21,DGM22,AMAT(1,iEQ))
        CALL DDMRADSUB2 (iEQ,C1,DGM21)      ! Constant term
        AMAT(jH   ,iEQ) = AMAT(jH   ,iEQ) + CX/Y(jH)
        AMAT(jSO4 ,iEQ) = AMAT(jSO4 ,iEQ) + CX/Y(jSO4)
        AMAT(jHSO4,iEQ) = AMAT(jHSO4,iEQ) - CX/Y(jHSO4)
CCC   ENDIF
C
C     Mass balance: NH3
C
      iEQ = iMBNH3
      IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XMAT(      iEQ,J) = SINI(jTNH3,J)
        ENDDO
        AMAT(jNH3aq ,iEQ) = 1. + 1./XLNH3
        AMAT(jNH4   ,iEQ) = 1.
      ENDIF
C
C     Mass balance: HNO3
C
      iEQ = iMBHNO3
      IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XMAT(      iEQ,J) = SINI(jTHNO3,J)
        ENDDO
        AMAT(jHNO3aq,iEQ) = 1. + 1./XLHNO3
        AMAT(jNO3   ,iEQ) = 1.
      ENDIF
C
C     Mass balance: FOA
C
      iEQ = iMBFOA
CCC   IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XMAT(      iEQ,J) = SINI(jTFOA,J)
        ENDDO
        AMAT(jFOAaq ,iEQ) = 1. + 1./XLFOA
        AMAT(jHCO2  ,iEQ) = 1.
CCC   ENDIF
C
C     Mass balance: CO2
C
      iEQ = iMBCO2
CCC   IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XMAT(      iEQ,J) = SINI(jTCO2,J) + SINI(jTCAC,J)
        ENDDO
        AMAT(jCO2aq ,iEQ) = 1. + 1./XLCO2
        AMAT(jHCO3  ,iEQ) = 1.
        AMAT(jCO3   ,iEQ) = 1.
CCC   ENDIF
C
C     Mass balance: SO2
C
      iEQ = iMBSO2
      IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XMAT(      iEQ,J) = SINI(jTSO2,J) ! Renewed
        ENDDO
        AMAT(jSO2aq ,iEQ) = 1. + 1./XLSO2
        AMAT(jHSO3  ,iEQ) = 1.
        AMAT(jSO3   ,iEQ) = 1.
      ENDIF
C
C     Mass balance: S(VI)
C
      iEQ = iMBSVI
CCC   IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XMAT(      iEQ,J) = SINI(jTSVI,J) ! Renewed
        ENDDO
        AMAT(jHSO4  ,iEQ) = 1.
        AMAT(jSO4   ,iEQ) = 1.
CCC   ENDIF
C
C     Charge Balance
C
      iEQ = iCB
CCC   IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XMAT(      iEQ,J) = -2.*SINI(jTCAC,J)
          ! Mg is assumed background constants
        ENDDO
        AMAT(jH   ,iEQ) =  1.
        AMAT(jNH4 ,iEQ) =  1.
        AMAT(jNO3 ,iEQ) = -1.
        AMAT(jHCO2,iEQ) = -1.
        AMAT(jHCO3,iEQ) = -1.
        AMAT(jHSO3,iEQ) = -1.
        AMAT(jHSO4,iEQ) = -1.
        AMAT(jOH  ,iEQ) = -1.
        AMAT(jCO3 ,iEQ) = -2.
        AMAT(jSO3 ,iEQ) = -2.
        AMAT(jSO4 ,iEQ) = -2.
CCC   ENDIF
C====================================================================
C
C     Solve the sensitivity eqns
C
      CALL DDMEQNSLV (NSENS,LSEQ,LSEN,AMAT,NDDMSP,XMAT,
     &                NDIM,AA,BB,IPVT,IPAR,IOUT)
      IF (IPAR.LT.0) THEN
        WRITE(IOUT,'(A)') ' Called by DDMRADEQL'
        RETURN
      ENDIF
C
C     Independent equilibria
C
      DO J = 1, NDDMSP
        IF (Y(jH2O2aq).LE.0.0) SINI(jTH2O2,J) = 0.0   ! If TH2O2 is depleted
        IF (Y(jO3aq  ).LE.0.0) SINI(jTO3  ,J) = 0.0   ! If TO3   is depleted
        IF (Y(jMHPaq ).LE.0.0) SINI(jTMHP ,J) = 0.0   ! If TMHP  is depleted
        IF (Y(jPAAaq ).LE.0.0) SINI(jTPAA ,J) = 0.0   ! If TPAA  is depleted
        IF (Y(jGLYaq ).LE.0.0) SINI(jTGLY ,J) = 0.0   ! If TGLY  is depleted
        IF (Y(jMGLYaq).LE.0.0) SINI(jTMGLY,J) = 0.0   ! If TMGLY is depleted
        IF (Y(jGLYDaq).LE.0.0) SINI(jTGLYD,J) = 0.0   ! If TGLYD is depleted
        ! RADM does not consume OH radical
        XMAT(jH2O2aq,J) = XLH2O2/(1.+XLH2O2) * SINI(jTH2O2,J) ! Renewed
        XMAT(jO3aq  ,J) = XLO3  /(1.+XLO3  ) * SINI(jTO3  ,J) ! Renewed
        XMAT(jMHPaq ,J) = XLMHP /(1.+XLMHP ) * SINI(jTMHP ,J) ! Renewed
        XMAT(jPAAaq ,J) = XLPAA /(1.+XLPAA ) * SINI(jTPAA ,J) ! Renewed
        XMAT(jGLYaq ,J) = XLGLY /(1.+XLGLY ) * SINI(jTGLY ,J) ! Renewed
        XMAT(jMGLYaq,J) = XLMGLY/(1.+XLMGLY) * SINI(jTMGLY,J) ! Renewed
        XMAT(jGLYDaq,J) = XLGLYD/(1.+XLGLYD) * SINI(jTGLYD,J) ! Renewed
        XMAT(jHOaq  ,J) = XLHO  /(1.+XLHO  ) * SINI(jTHO  ,J)
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMRADEQL


      SUBROUTINE DDMRADSUB1 (N0,N1,N2,C,DGM1,DGM2,VEC)
      INTEGER   N0,N1,N2,K
      REAL      C,DGM1,DGM2,VEC(N0)
      DO K = 1, N1
        VEC(K) = C * DGM1
      ENDDO
      DO K = N1 + 1, N2
        VEC(K) = C * DGM2
      ENDDO
      RETURN
      END
C     END OF SUBROUTINE DDMRADSUB1


      SUBROUTINE DDMRADSUB2 (I,C,DGM1)
      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmradm.inc'
      INTEGER   I,J
      REAL      C,DGM1
      DO J = 1, NDDMSP
        XMAT(I,J) = -C*DGM1*(
     &                       2.*SINI(jTNACL,J) ! (z1)^2 * (SINI(Na) + SINI(Cl))
     &                      +4.*SINI(jTCAC ,J) ! (z2)^2 *  SINI(Ca)
     &                       ) ! Fe, Mn, Mg & K are assumed background constants
      ENDDO
      RETURN
      END
C     END OF SUBROUTINE DDMRADSUB2


      SUBROUTINE DDMDAVIES (SION,GM1,GM2,DGM11,DGM12,DGM21,DGM22)

      REAL      LN10
      PARAMETER (LN10=2.3025851)

      REAL      SION, GM1, GM2, DGM11, DGM12, DGM21, DGM22
      REAL      DGMLOG, DELI1, DELI2

      DGMLOG = -0.509 * ( 0.5/(SION*(1.+SION)*(1.+SION)) - 0.2 )
      DELI1  = 0.5 ! dI/dY1 where Y1 = ions with absolute charge = 1
      DELI2  = 2.0 ! dI/dY2 where Y2 = ions with absolute charge = 2
      DGM11  = LN10 * GM1 * DGMLOG * DELI1    ! dGM1/dY1
      DGM12  = LN10 * GM1 * DGMLOG * DELI2    ! dGM1/dY2
      DGM21  = LN10 * GM2 * 4.*DGMLOG * DELI1 ! dGM2/dY1
      DGM22  = LN10 * GM2 * 4.*DGMLOG * DELI2 ! dGM2/dY2

      RETURN
      END
C     END OF SUBROUTINE DDMDAVIES


      SUBROUTINE DDMRADRXN (BB,SIV,TS6,MN,FE,
     &                      RH2O2,RO3A,RO3B,RO3C,RMHP,RPAA,
     &                      RGLY3,RMGLY3,RGLYD3,
     &                      YGLY,YMGLY,YGLYD,
     &                      DSIV_SCALE,DTW,IOUT)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmradm.inc'

      REAL      BB            ! Cloudwater pH
      REAL      SIV           ! [SO2(aq)] + [HSO3-] + [SO3=] (mol/liter)
      REAL      TS6           ! [SO4=] + [HSO4-]
      REAL      MN            ! Mn++ conc in cloudwater (mol/liter)
      REAL      FE            ! Fe+++ conc in cloudwater (mol/liter)
      REAL      RH2O2,RO3A,RO3B,RO3C,RMHP,RPAA
      REAL      RGLY3,RMGLY3,RGLYD3
      REAL      YGLY,YMGLY,YGLYD
      REAL      DSIV_SCALE
      REAL      DTW           ! Cloud chemistry timestep (sec)
      INTEGER   IOUT          ! Unit number of log file

      REAL      SH,C1,XSIV,XSVI
      REAL      R1XH,R1XH2O2,R1XHSO3
      REAL      R2XSO2,R2XHSO3,R2XSO3,R2XO3
      REAL      RO2A,RO2B,R3XSIV,R3XH,R3XSVI
      REAL      R4XH,R4XHSO3,R4XMHP
      REAL      R5XH,R5XHSO3,R5XPAA
      REAL      DR0,DR1,DR2,DR3,DR4,DR5,DR6,DR7,DR8 ! Sensitivities of rxn rates
      INTEGER   J
C
C     Initial check
C
      IF ( .NOT.LDDM ) RETURN             ! Check if DDM is being used
      IF ( IPAR.LE.0 ) RETURN             ! Check error from DDMRADEQL
C
C     Calculate coefficients
C
      SH      = Y(jH)*Y(jH)
C     R1
      C1      = 1.0 + 13.0*Y(jH)
      R1XH2O2 = Y(jHSO3  )*Y(jH)/C1
      R1XHSO3 = Y(jH2O2aq)*Y(jH)/C1
      R1XH    = Y(jH2O2aq)*Y(jHSO3  )/(C1*C1)
C     R2
      R2XSO2  = -RO3A*Y(jO3aq)
      R2XHSO3 = -RO3B*Y(jO3aq)
      R2XSO3  = -RO3C*Y(jO3aq)
      R2XO3   = -RO3A*Y(jSO2aq)-RO3B*Y(jHSO3)-RO3C*Y(jSO3)
C     R3
      RO2A    = -750.0*MN -2600.0*FE
      C1      = 1.0 + 75.0*(TS6**0.67)
      IF (BB.GE.4.2) THEN
        RO2B  = -2.51E13*MN*FE
        R3XSIV = (RO2A + RO2B*(Y(jH)**0.67))/C1
        R3XH   = (0.67*RO2B*SIV/(Y(jH)**0.33))/C1
      ELSE
        RO2B  = -3.72E7 *MN*FE
        R3XSIV = (RO2A + RO2B/(Y(jH)**0.74))/C1
        R3XH   = (-0.74*RO2B*SIV/(Y(jH)**1.74))/C1
      ENDIF
      R3XSVI = -R3XSIV*SIV*(0.67*75.0/(TS6**0.33))/C1
C     R4
      R4XH    = Y(jMHPaq)*Y(jHSO3)
      R4XHSO3 = Y(jH)*Y(jMHPaq)
      R4XMHP  = Y(jH)*Y(jHSO3)
C     R5
      C1      = -RPAA*Y(jH) -7.0E2
      R5XH    = -RPAA*Y(jPAAaq)*Y(jHSO3)
      R5XHSO3 = C1*Y(jPAAaq)
      R5XPAA  = C1*Y(jHSO3)

      DO J = 1, NDDMSP
        XSIV = XMAT(jSO2aq,J) + XMAT(jHSO3,J) + XMAT(jSO3,J)
        XSVI =                  XMAT(jHSO4,J) + XMAT(jSO4,J)
C
C     Oxidation by H2O2
C
        DR1 = -RH2O2 * ( R1XH2O2 * XMAT(jH2O2aq,J)
     &                 + R1XHSO3 * XMAT(jHSO3  ,J)
     &                 + R1XH    * XMAT(jH     ,J) )
C
C     Oxidation by O3
C
        DR2 = R2XSO2  * XMAT(jSO2aq,J)
     &      + R2XHSO3 * XMAT(jHSO3 ,J)
     &      + R2XSO3  * XMAT(jSO3  ,J)
     &      + R2XO3   * XMAT(jO3aq ,J)
C
C     Oxidation by O2 catalyzed by Fe and Mn ions
C
        DR3 = R3XSIV * XSIV
     &      + R3XH   * XMAT(jH  ,J)
     &      + R3XSVI * XSVI
C
C     Oxidation by MHP
C
        DR4 = -RMHP  * ( R4XH    * XMAT(jH     ,J)
     &                 + R4XHSO3 * XMAT(jHSO3  ,J)
     &                 + R4XMHP  * XMAT(jMHPaq ,J) )
C
C     Oxidation by PAA
C
        DR5 = R5XH    * XMAT(jH    ,J)
     &      + R5XHSO3 * XMAT(jHSO3 ,J)
     &      + R5XPAA  * XMAT(jPAAaq,J)
C
C     Oxidation of total S(IV)
C
        DR0 = DR1 + DR2 + DR3 + DR4 + DR5
C
C     GLY + OH
C
        DR6 = -RGLY3  * ( Y(jGLYaq ) * XMAT(jHOaq  ,J)
     &                  + Y(jHOaq  ) * XMAT(jGLYaq ,J) )
C
C     MGLY + OH
C
        DR7 = -RMGLY3 * ( Y(jMGLYaq) * XMAT(jHOaq  ,J)
     &                  + Y(jHOaq  ) * XMAT(jMGLYaq,J) )
C
C     GLYD + OH
C
        DR8 = -RGLYD3 * ( Y(jGLYDaq) * XMAT(jHOaq  ,J)
     &                  + Y(jHOaq  ) * XMAT(jGLYDaq,J) )
C
C     Update SINI
C
        SINI(jTH2O2,J) = SINI(jTH2O2,J) + DR1 * DTW * DSIV_SCALE
        SINI(jTO3  ,J) = SINI(jTO3  ,J) + DR2 * DTW * DSIV_SCALE
        SINI(jTMHP ,J) = SINI(jTMHP ,J) + DR4 * DTW * DSIV_SCALE
        SINI(jTPAA ,J) = SINI(jTPAA ,J) + DR5 * DTW * DSIV_SCALE
        SINI(jTSO2 ,J) = SINI(jTSO2 ,J) + DR0 * DTW * DSIV_SCALE
        SINI(jTSVI ,J) = SINI(jTSVI ,J) - DR0 * DTW * DSIV_SCALE
        SINI(jTGLY ,J) = SINI(jTGLY ,J) + DR6 * DTW
        SINI(jTMGLY,J) = SINI(jTMGLY,J) + DR7 * DTW
        SINI(jTGLYD,J) = SINI(jTGLYD,J) + DR8 * DTW
        SINI(jTSOAC,J) = SINI(jTSOAC,J) - (YGLY*DR6 + YMGLY*DR7 + YGLYD*DR8) * DTW
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMRADRXN

