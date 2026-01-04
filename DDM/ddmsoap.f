C=====================================================================
C
C   DDM-PM SUBROUTINES FOR SOAP
C
C   Subroutines included in this file:
C     DDMSOAINI
C     DDMSOAP
C     DDMSOATRV (trivial solutions)
C     DDMSOAAP (absorptive partitioning)
C
C   Dependent upon:
C     tracer.com -> camx.prm  - for MXFDDM, NDDMSP & LDDM
C     ddmsoap.inc
C     soap.inc
C
C   Log:
C     05/20/2005  bkoo - original development
C     06/25/2008  bkoo - revised for CAMx v4.5 SOA module update
C     03/18/2014  bkoo - revised for benzene SOA
C     08/25/2016  bkoo - updated for new SOAP
C     10/09/2018  gy   - updated for aerosol photolysis
C
C=====================================================================

      SUBROUTINE DDMSOAINI (CONVFAC,SDDM,NFAM,NSEN)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmsoap.inc'
      INCLUDE 'soap.inc'

      INTEGER   NFAM, NSEN
      REAL      SDDM(NFAM,NSEN)           ! [ppm]/[p] for gas; [ug/m3]/[p] for PM
      REAL      CONVFAC                   ! umol/m3 = ppm * convfac

      INTEGER   I, J
      logical, save :: lfirsttime = .true.

      if (lfirsttime) then
        lfirsttime = .false.
        if (NSENS.ne.NSOAP) stop'ERROR in DDMSOAINI: NSENS .ne. NSOAP'
      endif

      DO J = 1, NFAM
        SIPA(J) = SDDM(J,jdpSOPA)
        SIPB(J) = SDDM(J,jdpSOPB)
        SPRE(J) = SDDM(J,jdpPOA )/MWPOA
     &          + SDDM(J,jdpSOPA)/MWSOPA
     &          + SDDM(J,jdpSOPB)/MWSOPB  ! [mole]/[p]
        DO I = 1, NSENS
          SINI(I,J) = SDDM(J,jdpCG1 +I-1)*CONVFAC*MWSOAP(I)
     &              + SDDM(J,jdpSOA1+I-1) ! [ug/m3]/[p]
          XMAT(I,J) = 0.0                 ! Initialize the solution matrix
        ENDDO
      ENDDO

      IPAR = 1

      RETURN
      END
C     END OF SUBROUTINE DDMSOAINI


      SUBROUTINE DDMSOAP (DT,AJNO2,CONVFAC,SDDM,NFAM,NSEN,IERR,IOUT)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmsoap.inc'
      INCLUDE 'soap.inc'

      INTEGER   NFAM, NSEN
      REAL      SDDM(NFAM,NSEN)           ! [ppm]/[p] for gas; [ug/m3]/[p] for PM
      REAL      CONVFAC                   ! umol/m3 = ppm * convfac
      REAL      DT                        ! Time step [hours]
      REAL      AJNO2                     ! Accumulated NO2 photolysis rate ([1/hr]*[hr])
      INTEGER   IERR, IOUT

      REAL      XGAS(NSENS,MXFDDM)        ! Solution matrix for gases
      REAL      FPOLYA, FPOLYB, SPOLY, FPHOTR
      INTEGER   I, J
C
C     Initialize error flag
C
      IERR = 1
C
C     Check error from the previous routines
C
      IF ( IPAR.LT.0 ) THEN
        IERR = IPAR
        RETURN
      ENDIF
C
C     Calculate the sensitivities of gas species
C
      DO J = 1, NDDMSP
        DO I = 1, NSENS
          XGAS(I,J) = SINI(I,J) - XMAT(I,J)
        ENDDO
      ENDDO
C
C     Calculate changes due to aerosol photolysis followed by polymerization
C
      FPHOTR = EXP( -SFACJSOA * AJNO2 )   ! Fraction of remaining SOA after photolysis
      FPOLYA = 1.0 - EXP( -KPOLYA * DT )  ! polymerized fraction for anthro SOA
      FPOLYB = 1.0 - EXP( -KPOLYB * DT )  ! polymerized fraction for bio SOA
      DO J = 1, NDDMSP
        SIPA(J) = FPHOTR * SIPA(J)
        SIPB(J) = FPHOTR * SIPB(J)
        DO I = 1, NSOAA
          XMAT(I,J) = FPHOTR * XMAT(I,J)
          SPOLY = FPOLYA * XMAT(I,J)
          XMAT(I,J) = XMAT(I,J) - SPOLY
          SIPA(J) = SIPA(J) + SPOLY
        ENDDO
        DO I = NSOAA+1, NSENS
          XMAT(I,J) = FPHOTR * XMAT(I,J)
          SPOLY = FPOLYB * XMAT(I,J)
          XMAT(I,J) = XMAT(I,J) - SPOLY
          SIPB(J) = SIPB(J) + SPOLY
        ENDDO
      ENDDO
C
C     Assign the solution
C
      DO J = 1, NFAM
        DO I = 1, NSENS
          SDDM(J,jdpCG1 +I-1) = XGAS(I,J)/MWSOAP(I)/CONVFAC
          SDDM(J,jdpSOA1+I-1) = XMAT(I,J)
        ENDDO
        SDDM(J,jdpSOPA) = SIPA(J)
        SDDM(J,jdpSOPB) = SIPB(J)
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMSOAP


      SUBROUTINE DDMSOATRV (CTOT,CSAT,CPRE,MSOA,CAER,IDX,ICASE,IOUT)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmsoap.inc'

      REAL      CTOT
      REAL      CSAT
      REAL      CPRE          ! Only for CASE 3 (pre-existing OA in mole)
      REAL      MSOA          ! Only for CASE 3 (MW_soa)
      REAL      CAER          ! Only for CASE 3
      INTEGER   IDX           ! Index in the base set of SOA species
      INTEGER   ICASE         ! Case ID: must be 1, 2, or 3
      INTEGER   IOUT          ! Unit number of log file

      REAL      C1, C2, C3
      INTEGER   J
C     
C     Initial check
C     
      IF ( .NOT.LDDM ) RETURN             ! Check if DDM is being used

      GOTO (100, 200, 300) ICASE
      WRITE(IOUT,'(//,A)') 'ERROR in DDMSOATRV:'
      WRITE(IOUT,'(A,I5)')  ' Invalid ICASE - ',ICASE
      IPAR = -1
      RETURN
C
C     CASE 1: non-solution-forming case OR pure component w/o pre-existing OA
C
100   CONTINUE
      IF (CTOT.LE.CSAT) GOTO 200
      DO J = 1, NDDMSP
        XMAT(IDX,J) = SINI(IDX,J)
      ENDDO
      RETURN
C
C     CASE 2: insignificant total concentration OR below threshold for partitioning
C
200   CONTINUE
      DO J = 1, NDDMSP
        XMAT(IDX,J) = 0.0
      ENDDO
      RETURN
C
C     CASE 3: analytical solution of a trivial case
C
300   CONTINUE
      C1 = CAER + MSOA * CPRE
      C2 = CTOT - CAER
      C3 = C1 + CSAT - C2
      IF (C3.EQ.0.0) THEN
        WRITE(IOUT,'(//,A)') 'ERROR in DDMSOATRV:'
        WRITE(IOUT,'(A)')    ' Zero denominator!'
        IPAR = -1
        RETURN
      ENDIF
      DO J = 1, NDDMSP
        XMAT(IDX,J) = ( C1*SINI(IDX,J) + MSOA*C2*SPRE(J) ) / C3
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMSOATRV


      SUBROUTINE DDMSOAAP (CTOT,CSAT,MSOA,TOM,NSOL,IDX,IOUT)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmsoap.inc'

      REAL      CTOT(NSENS)
      REAL      CSAT(NSENS)
      REAL      MSOA(NSENS)
      REAL      TOM           ! Total Organic Matter (moles)
      INTEGER   NSOL, IDX(NSENS)
      INTEGER   IOUT          ! Unit number of log file

      REAL      DENOM, DTOM
      REAL      SUM1, SUM2, T1(NSENS), T2(NSENS)
      INTEGER   I, J
C     
C     Initial check
C     
      IF ( .NOT.LDDM ) RETURN             ! Check if DDM is being used
C
C     Calculate the sensitivity of TOM
C
      DO J = 1, NDDMSP
        SUM1 = 0.0
        SUM2 = 0.0
        DO I = 1, NSOL
          DENOM = MSOA(I)*TOM + CSAT(I)
          T1(I) = SINI(IDX(I),J) * TOM / DENOM
          T2(I) = CTOT(I)*CSAT(I) / (DENOM*DENOM)
          SUM1  = SUM1 + T1(I)
          SUM2  = SUM2 + T2(I)
        ENDDO
        IF (SUM2.EQ.1.0) THEN
          WRITE(IOUT,'(//,A)') 'ERROR in DDMSOAAP:'
          WRITE(IOUT,'(A)')    ' Zero denominator!'
          IPAR = -1
          RETURN
        ENDIF
        DTOM = (SUM1 + SPRE(J)) / (1. - SUM2)
        DO I = 1, NSOL
          XMAT(IDX(I),J) = (T1(I) + T2(I)*DTOM) * MSOA(I)
        ENDDO
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMSOAAP

