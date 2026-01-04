C=====================================================================
C
C   MATRIX SOLVER SUBROUTINES FOR DDM-PM
C
C   Subroutines included in this file:
C     DDMEQNSUB
C     DDMEQNSLV
C
C   Dependent upon:
C     none
C
C   Log:
C     03/03/2005  bkoo - original development
C
C=====================================================================

      SUBROUTINE DDMEQNSUB (L1,L2)
      LOGICAL   L1,L2
      L1 = .FALSE.
      L2 = .FALSE.
      RETURN
      END
C     END OF SUBROUTINE DDMEQNSUB


      SUBROUTINE DDMEQNSLV (NSENS,LSEQ,LSEN,AMAT,NDDMSP,XLIQ,
     &                      NDIM,AA,BB,IPVT,IERR,IOUT)

      INTEGER   NSENS, NDDMSP, NDIM
      LOGICAL   LSEQ(NSENS)               ! Logical flag for eqns
      LOGICAL   LSEN(NSENS)               ! Logical flag for sensitivities
      REAL      AMAT(NSENS,NSENS),XLIQ(NSENS,NDDMSP)

      REAL      AA(NDIM,NDIM),BB(NDIM)    ! Now adjustable arrays
      INTEGER   IPVT(NDIM)
      INTEGER   IERR, IOUT

      INTEGER   IDX(NDIM)
      INTEGER   IEQN, JEQN, IVAR, JVAR, J, INFO, JOB
      INTEGER   NERR, MXERR
      PARAMETER( MXERR = 5 )              ! Max error count
      LOGICAL   LAFAC
      REAL      AFAC, RLMT
      PARAMETER( AFAC = 1.E-2 )           ! Adjustment factor
C
C     Initialization
C
      NERR = 0
      RLMT = 1.E37                        ! Limiting magnitude of the coeffs
C
C     Fill AA and BB
C
      JEQN = 0
      DO IEQN = 1, NSENS
        IF ( LSEQ(IEQN) ) THEN
          JEQN = JEQN + 1
          JVAR = 0
          DO IVAR = 1, NSENS
            IF ( LSEN(IVAR) ) THEN
              JVAR = JVAR + 1
              AA(JVAR,JEQN) = AMAT(IVAR,IEQN)
            ENDIF
          ENDDO
          IDX(JEQN) = IEQN                ! For BB
        ENDIF
      ENDDO
C
C     Adjust the magnitude of the coefficients to avoid "Floating exception" error
C
100   CONTINUE
      DO JEQN = 1, NDIM
        LAFAC = .FALSE.
        DO JVAR = 1, NDIM
          IF ( ABS(AA(JVAR,JEQN)) .GT. RLMT ) THEN
            LAFAC = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF ( LAFAC ) THEN
          DO JVAR = 1, NDIM
            AA(JVAR,JEQN) = AA(JVAR,JEQN) * AFAC
          ENDDO
          DO J = 1, NDDMSP
            XLIQ(IDX(JEQN),J) = XLIQ(IDX(JEQN),J) * AFAC
          ENDDO
        ENDIF
      ENDDO
C
C     Decompose matrix AA - check for zero determinant
C
      CALL SGEFA (AA,NDIM,NDIM,IPVT,INFO)
      IF (INFO.NE.0) THEN
        NERR = NERR + 1
        RLMT = 1.E11
        IF (NERR.LE.MXERR) GOTO 100
        WRITE(IOUT,'(//,A)') 'ERROR in DDMEQNSLV:'
        WRITE(IOUT,'(A)')    ' Zero determinant in SGEFA!'
        WRITE(IOUT,*)        ' info = ',INFO
        IERR = -1
        RETURN
      ENDIF
C
C     Solve the system
C
      JOB = 1                             ! Solve trans(AA) * X = BB
      DO J = 1, NDDMSP
        DO JEQN = 1, NDIM
          BB(JEQN) = XLIQ(IDX(JEQN),J)
        ENDDO
        CALL SGESL (AA,NDIM,NDIM,IPVT,BB,JOB)
C
C     Map the answers
C
        JVAR = 0
        DO IVAR = 1, NSENS
          XLIQ(IVAR,J) = 0.0              ! Clearing
          IF ( LSEN(IVAR) ) THEN
            JVAR = JVAR + 1
            XLIQ(IVAR,J) = BB(JVAR)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMEQNSLV

