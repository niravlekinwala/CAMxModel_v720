C=====================================================================
C
C   DDM-PM SUBROUTINES FOR ISORROPIA
C
C   Subroutines included in this file:
C     DDMISRINI
C     DDMISRPIA
C     DDMISRDRY
C     DDMISRLIQ
C     DDMDIFACT
C     DDMKMFUL
C     DDMMKBI
C     DDMDIFZSR
C     DDMISRMDR
C     DDMCALCHA
C     DDMCALCNA
C     DDMCALCNHA
C     DDMCALCNH3
C     DDMISRSAVH
C
C   Dependent upon:
C     tracer.com -> camx.prm  - for MXFDDM, NDDMSP & LDDM
C     ddmisrpia.inc
C     isrpia.inc
C
C   Log:
C     03/03/2005  bkoo - original development
C     06/25/2008  bkoo - changes due to ISORROPIA v1.7 (03/26/07):
C                        1. LN10 -> double precision (DDMDIFACT)
C                        2. upper bound of IONIC changed to 500 (DDMDIFACT) -> 100
C                        3. GAMA bounded between 10^-11 and 10^11 (DDMDIFACT) -> 10^-5 ~ 10^5
C                        4. temperature correction bug fix (DDMKMFUL)
C     09/19/2013  bkoo - revised for excluding mineral dust nitrate
C     12/07/2014  bkoo - several modifications to enhance robustness
C     10/26/2018  bkoo - updated for revised Ca scaling
C
C=====================================================================

      SUBROUTINE DDMISRINI (LNO3ISR,fCA_PRM,fCA_CRS,
     &                      CONVFAC,WTMOL,SDDM,NFAM,NSEN,NWTMOL,IOUT)

      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      INTEGER   NFAM, NSEN, NWTMOL, IOUT
      REAL      SDDM(NFAM,NSEN)           ! [ppm]/[p] for gas; [ug/m3]/[p] for PM
      REAL      WTMOL(NWTMOL)             ! Molecular weight array
      REAL      fCA_PRM,fCA_CRS           ! mass fraction of Ca on FPRM/FCRS
      REAL      CONVFAC                   ! umol/m3 = ppm * convfac
      LOGICAL   LNO3ISR                   ! flag to indicate non-zero nitrate for isorropia

      INTEGER   I, J

      DO J = 1, NFAM
C     LOADDMPM has initialized SDDM rows for non-model species to zero
        SINI(jTNA ,J) =   SDDM(J,jdpNA  )/WTMOL(5)  *1.E-6
        SINI(jTSO4,J) = ( SDDM(J,jdpSULF)*CONVFAC + 
     &                    SDDM(J,jdpPSO4)/WTMOL(1) )*1.E-6
        SINI(jTNH4,J) = ( SDDM(J,jdpNH3 )*CONVFAC + 
     &                    SDDM(J,jdpPNH4)/WTMOL(2) )*1.E-6
        SINI(jTNO3,J) = ( SDDM(J,jdpHNO3)*CONVFAC + 
     &                    SDDM(J,jdpPNO3)/WTMOL(3) )*1.E-6
        IF ( LNO3ISR ) THEN
          SINI(jTNO3,J) = SINI(jTNO3,J) -
     &                  ( fCA_PRM*SDDM(J,jdpFPRM)
     &                  + fCA_CRS*SDDM(J,jdpFCRS) )/WTMOL(8)*1.E-6
        ELSE
          SINI(jTNO3,J) = 0.0
        ENDIF
        SINI(jTCL ,J) = ( SDDM(J,jdpHCL )*CONVFAC + 
     &                    SDDM(J,jdpPCL )/WTMOL(4) )*1.E-6
        DO I = 1, NSENS
          XMAT(I,J) = 0.0                 ! Initialize the solution matrix
        ENDDO
      ENDDO
      LMDRH = .FALSE.                     ! Initialize LMDRH flag
      IPAR = IOUT                         ! Unit number of log file

      RETURN
      END
C     END OF SUBROUTINE DDMISRINI


      SUBROUTINE DDMISRPIA (LNO3ISR,fCA_PRM,fCA_CRS,CONVFAC,
     &                      WTMOL,SDDM,NFAM,NSEN,NWTMOL,IERR,IOUT)

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      INTEGER   NFAM, NSEN, NWTMOL, IERR, IOUT
      REAL      SDDM(NFAM,NSEN)           ! [ppm]/[p] for gas; [ug/m3]/[p] for PM
      REAL      WTMOL(NWTMOL)             ! Molecular weight array
      REAL      fCA_PRM,fCA_CRS           ! mass fraction of Ca on FPRM/FCRS
      REAL      CONVFAC                   ! umol/m3 = ppm * convfac
      LOGICAL   LNO3ISR                   ! flag to indicate non-zero nitrate for isorropia

CCC   REAL*8    TMPH
      CHARACTER SC*1                      ! The first character of SCASE
      INTEGER   J
C
C     Initialize error flag
C
      IERR = 1
C
C     Check for a trivial solution
C
      SC = SCASE(1:1)
      IF (SC.EQ.'?') THEN                 ! No aerosol mass
        IERR = 0                          ! -> no change in sensitivities
        RETURN
      ENDIF
C
C     MDRH region
C
      IF ( LMDRH ) THEN
        IF ( IPAR.LT.0 ) THEN             ! Check error from DDMISRMDR
          IERR = IPAR
          RETURN
        ENDIF
CCC     IF (SC.EQ.'B') THEN
CCC       CALL DDMCALCNH3
CCC     ELSEIF (SC.EQ.'E') THEN
CCC       TMPH      = MOLAL(jH)
CCC       MOLAL(jH) = SAVH
CCC       CALL DDMCALCNA
CCC       MOLAL(jH) = TMPH
CCC       CALL DDMCALCNH3
CCC     ELSEIF (SC.EQ.'I') THEN
CCC       TMPH      = MOLAL(jH)
CCC       MOLAL(jH) = SAVH
CCC       CALL DDMCALCNHA
CCC       MOLAL(jH) = TMPH
CCC       CALL DDMCALCNH3
CCC     ENDIF
CCC   This is not always working - Let's make it simple
        GOTO 100
      ENDIF
C
C     Solid only phase
C
      IF ( WATER.LE.TINY ) THEN           ! Dry

10      CONTINUE
        DO J = 1, NDDMSP
          CALL DDMISRDRY (SINI(1,J),XMAT(1,J),IERR,IOUT)
          IF (IERR.LT.0) RETURN
        ENDDO

        GOTO 100
      ENDIF
C
C     Liquid & solid phase
C
      IF ( IONIC.LT.TINY2 ) GOTO 10

CCC      MOLAL(jH     ) = MAX(MOLAL(jH     ), TINY)
CCC      MOLAL(jHSO4  ) = MAX(MOLAL(jHSO4  ), TINY)
      CALL DDMISRLIQ (XMAT,IERR,IOUT)
      IF (IERR.LT.0) RETURN

100   CONTINUE
C
C     Assign the solution
C
      DO J = 1, NFAM
        SDDM(J,jdpSULF) =  0.0
        SDDM(J,jdpPSO4) =  SINI(jTSO4,J)*WTMOL(1)*1.E6
        SDDM(J,jdpNH3 ) =  XMAT(jNH3 ,J)/CONVFAC*1.E6
        SDDM(J,jdpPNH4) = (SINI(jTNH4,J) - XMAT(jNH3 ,J))*WTMOL(2)*1.E6
        IF ( LNO3ISR ) THEN
          SDDM(J,jdpHNO3) =  XMAT(jHNO3,J)/CONVFAC*1.E6
          SDDM(J,jdpPNO3) = (SINI(jTNO3,J) - XMAT(jHNO3,J)
     &                       + (fCA_PRM*SDDM(J,jdpFPRM)
     &                         +fCA_CRS*SDDM(J,jdpFCRS))/WTMOL(8)*1.E-6)
     &                                                  *WTMOL(3)*1.E6
        ELSE
          SDDM(J,jdpPNO3) = SDDM(J,jdpPNO3)
     &                    + SDDM(J,jdpHNO3)*CONVFAC*WTMOL(3)
          SDDM(J,jdpHNO3) = 0.0
        ENDIF
        SDDM(J,jdpHCL ) =  XMAT(jHCL ,J)/CONVFAC*1.E6
        SDDM(J,jdpPCL ) = (SINI(jTCL ,J) - XMAT(jHCL ,J))*WTMOL(4)*1.E6
        SDDM(J,jdpPH2O) =  XMAT(jH2O ,J)*1.E9
        SDDM(J,jdpPH2O) = MIN( MAX( SDDM(J,jdpPH2O), -1.E3 ), 1.E3 ) ! cbk: cap water sens
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMISRPIA


      SUBROUTINE DDMISRDRY (SI,XVEC,IERR,IOUT)

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL      SI(NINPT)                 ! Initial Sensitivities of total
      REAL      XVEC(NSENS)               ! Solution vector X
      INTEGER   IERR, IOUT

      CHARACTER SC*1
      REAL*8    OMPS, ALF, THETA2, FRNA, FRNO3, FRCL, FRSO4
      REAL      DOMPS, DZE, DALF, DFRNA, DFRNO3, DFRCL

      SC = SCASE(1:1)
      IF     (SC.EQ.'A') THEN
        XVEC(jNH42S4) = SI(jTSO4)
        XVEC(jNH3)    = SI(jTNH4) -2.*SI(jTSO4)
      ELSEIF (SC.EQ.'B' .OR. SC.EQ.'E') THEN
        IF (3.D0*W(jTSO4) .LE. 2.D0*W(jTNH4)) THEN
          XVEC(jLC)     = 2.*SI(jTSO4) -   SI(jTNH4)
          XVEC(jNH42S4) = 2.*SI(jTNH4) -3.*SI(jTSO4)
        ELSE
          XVEC(jLC)     =    SI(jTNH4) -   SI(jTSO4)
          XVEC(jNH4HS4) = 3.*SI(jTSO4) -2.*SI(jTNH4)
        ENDIF
        IF (SC.EQ.'E') XVEC(jHNO3) = SI(jTNO3)
      ELSEIF (SC.EQ.'D') THEN
        XVEC(jNH42S4) = SI(jTSO4)
        OMPS = W(jTNH4) -2.D0*W(jTSO4) -W(jTNO3)
        DOMPS=SI(jTNH4) -2. *SI(jTSO4)-SI(jTNO3)
        IF (OMPS.GE.ZERO) THEN
          DZE = -DOMPS*GHNO3 / (2.D0*GHNO3 +OMPS)
          IF (CNH4NO3.LE.ZERO) DZE = SI(jTNO3)
          XVEC(jNH4NO3) = SI(jTNO3) -DZE
          XVEC(jNH3)    = DOMPS     +DZE
          XVEC(jHNO3)   =            DZE
        ELSE
          DZE = DOMPS*GNH3 / (2.D0*GNH3 -OMPS)
          IF (CNH4NO3.LE.ZERO) DZE = SI(jTNH4) -2.*SI(jTSO4)
          XVEC(jNH4NO3) = SI(jTNH4) -2.*SI(jTSO4) -DZE
          XVEC(jNH3)    =                          DZE
          XVEC(jHNO3)   = -DOMPS                  +DZE
        ENDIF
      ELSEIF (SC.EQ.'G') THEN
        XVEC(jNA2SO4) = .5*SI(jTNA)
        XVEC(jNH42S4) = SI(jTSO4) -XVEC(jNA2SO4)
        ALF = W(jTNH4) -2.D0*CNH42S4
        DALF=SI(jTNH4) -2.*XVEC(jNH42S4)
        IF (CNH4NO3.GT.ZERO .AND. CNH4CL.GT.ZERO) THEN
          THETA2 = XK10/XK6
          XVEC(jNH4CL)  = SI(jTCL) + ( DALF-SI(jTNO3)-SI(jTCL) )*GHCL /
     &                  ( ALF-W(jTNO3)-W(jTCL)+2.D0*(ONE+THETA2)*GHCL )
          XVEC(jNH4NO3) = SI(jTNO3) - THETA2*( SI(jTCL)-XVEC(jNH4CL) )
        ELSEIF (CNH4CL.GT.ZERO) THEN
          XVEC(jNH4CL)  = (DALF*GHCL + SI(jTCL)*GNH3)/(GHCL + GNH3)
          XVEC(jNH4NO3) = 0.0
        ELSEIF (CNH4NO3.GT.ZERO) THEN
          XVEC(jNH4CL)  = 0.0
          XVEC(jNH4NO3) = (DALF*GHNO3 + SI(jTNO3)*GNH3)/(GHNO3 + GNH3)
        ELSE
          XVEC(jNH4CL)  = 0.0
          XVEC(jNH4NO3) = 0.0
        ENDIF
        XVEC(jNH3)    = DALF     -XVEC(jNH4CL) -XVEC(jNH4NO3)
        XVEC(jHCL)    = SI(jTCL) -XVEC(jNH4CL)
        XVEC(jHNO3)   = SI(jTNO3)              -XVEC(jNH4NO3)
      ELSEIF (SC.EQ.'H') THEN
        XVEC(jNA2SO4) = SI(jTSO4)
        FRNA = W(jTNA) -2.D0*W(jTSO4)
        DFRNA=SI(jTNA) -2. *SI(jTSO4)
        IF (FRNA.GE.W(jTNO3)) THEN
          XVEC(jNANO3)  = SI(jTNO3)
          FRNO3 = ZERO
          DFRNO3= 0.0
          IF (FRNA-W(jTNO3).GE.W(jTCL)) THEN
            XVEC(jNACL) = SI(jTCL)
            FRCL = ZERO
            DFRCL= 0.0
          ELSE
            XVEC(jNACL) = DFRNA -SI(jTNO3)
            FRCL = W(jTCL) +W(jTNO3) -FRNA
            DFRCL=SI(jTCL)+SI(jTNO3)-DFRNA
          ENDIF
        ELSE
          XVEC(jNANO3)  = DFRNA
          FRNO3 = W(jTNO3) -FRNA
          DFRNO3=SI(jTNO3)-DFRNA
          FRCL  = W(jTCL)
          DFRCL =SI(jTCL)
        ENDIF
        IF (CNH4NO3.GT.ZERO .AND. CNH4CL.GT.ZERO) THEN
          THETA2 = XK10/XK6
          XVEC(jNH4CL)  = DFRCL + ( SI(jTNH4)-DFRNO3-DFRCL )*GHCL /
     &                  ( W(jTNH4)-FRNO3-FRCL+2.D0*(ONE+THETA2)*GHCL )
          XVEC(jNH4NO3) = DFRNO3 - THETA2*( DFRCL-XVEC(jNH4CL) )
        ELSEIF (CNH4CL.GT.ZERO) THEN
          XVEC(jNH4CL)  = (SI(jTNH4)*GHCL + DFRCL*GNH3)/(GHCL + GNH3)
          XVEC(jNH4NO3) = 0.0
        ELSEIF (CNH4NO3.GT.ZERO) THEN
          XVEC(jNH4CL)  = 0.0
          XVEC(jNH4NO3) = (SI(jTNH4)*GHNO3 + DFRNO3*GNH3)/(GHNO3+GNH3)
        ELSE
          XVEC(jNH4CL)  = 0.0
          XVEC(jNH4NO3) = 0.0
        ENDIF
        XVEC(jNH3)    = SI(jTNH4) -XVEC(jNH4CL) -XVEC(jNH4NO3)
        XVEC(jHCL)    = DFRCL     -XVEC(jNH4CL)
        XVEC(jHNO3)   = DFRNO3                  -XVEC(jNH4NO3)
      ELSEIF (SC.EQ.'I') THEN
        XVEC(jNA2SO4) = .5*SI(jTNA)
        FRSO4 = W(jTSO4) -0.5D0*W(jTNA) -2.D0*W(jTNH4)/3.D0
        IF (FRSO4.LE.TINY) THEN
          XVEC(jLC)     = 2.*SI(jTSO4) -SI(jTNA) -SI(jTNH4)
          XVEC(jNH42S4) = 2.*SI(jTNH4) -3.*SI(jTSO4) +1.5*SI(jTNA)
        ELSE
          IF (FRSO4.LE.W(jTNH4)/3.D0) THEN
            XVEC(jNH4HS4) = 3.*SI(jTSO4) -1.5*SI(jTNA) -2.*SI(jTNH4)
            XVEC(jLC)     = SI(jTNH4) -SI(jTSO4) +.5*SI(jTNA)
          ELSE
            XVEC(jNH4HS4) = SI(jTNH4)
            IF (0.5D0*W(jTNA).GT.TINY) THEN
              XVEC(jNAHSO4) = 2.*SI(jTSO4) -SI(jTNA) -2.*SI(jTNH4)
              XVEC(jNA2SO4) = SI(jTNH4) +SI(jTNA) -SI(jTSO4)
            ENDIF
          ENDIF
        ENDIF
        XVEC(jHNO3) = SI(jTNO3)
        XVEC(jHCL)  = SI(jTCL)
      ELSE
        WRITE(IOUT,'(//,A)') 'ERROR in DDMISRDRY:'
        WRITE(IOUT,'(A,A)')  ' Case not supported - ',SCASE
        IERR = -1
        RETURN
      ENDIF

      RETURN
      END
C     END OF SUBROUTINE DDMISRDRY


      SUBROUTINE DDMISRLIQ (XLIQ,IERR,IOUT)

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL      XLIQ(NSENS,MXFDDM)        ! Solution matrix for wet aerosols
      INTEGER   IERR, IOUT

      CHARACTER SC*1
      LOGICAL   LSEQ(NSENS)               ! If true, the eqn is solved
      LOGICAL   LSEN(NSENS)               ! If true, the sensitivity is solved
      REAL      DELGAMA(NIONSP1,NPAIR)    ! del GAMA by del Y(A)
      REAL      AMAT(NSENS,NSENS)         ! Left-hand side matrix A
                                          ! Right-hand side vector B -> XLIQ
      REAL      AA(NSENS,NSENS),BB(NSENS) ! These are passed to DDMEQNSLV solely
      INTEGER   IPVT(NSENS)               ! to allow their size to be adjustable
      LOGICAL   LALFA                     ! If true, calculate DALFA
      REAL*8    CX, C1, C2
CCC   REAL*8    A1, A2, A3, XSUM1, XSUM2, XSUM3
      INTEGER   I, J, K, NDIM, iEQ

      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A__1, A__2, A__3, A__4, A__5, A__6, A__7, A__8
c$omp threadprivate(/SOLUT/)
C
C     Initialize arrays
C
      DO I = 1, NSENS
        LSEQ(I) = .TRUE.
        LSEN(I) = .TRUE.
        DO J = 1, NSENS
          AMAT(J,I) = 0.0
        ENDDO
      ENDDO
      LALFA = .FALSE.
C
C     Select eqns & sensitivities to solve
C
      SC = SCASE(1:1)
      IF     (SC.EQ.'A' .OR. SC.EQ.'B' .OR. SC.EQ.'C') THEN ! Na,Cl,NO3=0
        CALL DDMEQNSUB (LSEQ(iMBNA ),LSEN(jNA  ))
        CALL DDMEQNSUB (LSEQ(iMBCL ),LSEN(jCL  ))
        CALL DDMEQNSUB (LSEQ(iMBNO3),LSEN(jNO3 ))
        CALL DDMEQNSUB (LSEQ(iK3   ),LSEN(jHCL ))
        CALL DDMEQNSUB (LSEQ(iK4   ),LSEN(jHNO3))
      ELSEIF (SC.EQ.'D' .OR. SC.EQ.'E' .OR. SC.EQ.'F') THEN ! Na,Cl=0
        CALL DDMEQNSUB (LSEQ(iMBNA ),LSEN(jNA  ))
        CALL DDMEQNSUB (LSEQ(iMBCL ),LSEN(jCL  ))
        CALL DDMEQNSUB (LSEQ(iK3   ),LSEN(jHCL ))
      ELSEIF (SC.EQ.'G') THEN
        IF (SCASE.EQ.'G1 ; SUBCASE 2' .OR.
     &      SCASE.EQ.'G2 ; SUBCASE 1' .OR.
     &      SCASE.EQ.'G2 ; SUBCASE 3' .OR.
     &      SCASE.EQ.'G3 ; SUBCASE 1') THEN
          CALL DDMEQNSUB (LSEQ(iK5 ),LSEN(jNA  ))
          CNA2SO4 = 0.5D0*W(jTNA)
          MOLAL(jSO4) = MAX(MOLAL(jSO4) - 0.5D0*MOLAL(jNA ), ZERO)
          MOLAL(jNA ) = ZERO
        ENDIF
      ELSEIF (SC.EQ.'H') THEN
        IF (PSI5.LT.2.D0*TINY .OR. PSI4.LT.5.D-9)
     &                     CALL DDMEQNSUB (LSEQ(iK2   ),LSEN(jNH4 ))
        IF (GHNO3.LE.TINY) CALL DDMEQNSUB (LSEQ(iK4   ),LSEN(jHNO3)) ! Na takes all HNO3?
CCC      ELSEIF (SC.EQ.'I' .OR. SC.EQ.'J') THEN                ! CALCNHA can give TINY to NO3-/CL-
CCC        IF (MOLAL(jNO3).LE.TINY) CALL DDMEQNSUB (LSEQ(iK4),LSEN(jNO3))
CCC        IF (MOLAL(jCL ).LE.TINY) CALL DDMEQNSUB (LSEQ(iK3),LSEN(jCL ))
      ENDIF
      IF (CNH42S4.LE.ZERO) CALL DDMEQNSUB (LSEQ(iK7 ),LSEN(jNH42S4))
      IF (CNH4HS4.LE.ZERO) CALL DDMEQNSUB (LSEQ(iK12),LSEN(jNH4HS4))
      IF (CNACL  .LE.ZERO) CALL DDMEQNSUB (LSEQ(iK8 ),LSEN(jNACL  ))
      IF (CNA2SO4.LE.ZERO) CALL DDMEQNSUB (LSEQ(iK5 ),LSEN(jNA2SO4))
      IF (CNANO3 .LE.ZERO) CALL DDMEQNSUB (LSEQ(iK9 ),LSEN(jNANO3 ))
      IF (CNH4NO3.LE.ZERO) CALL DDMEQNSUB (LSEQ(iK10),LSEN(jNH4NO3))
      IF (CNH4CL .LE.ZERO) CALL DDMEQNSUB (LSEQ(iK6 ),LSEN(jNH4CL ))
      IF (CNAHSO4.LE.ZERO) CALL DDMEQNSUB (LSEQ(iK11),LSEN(jNAHSO4))
      IF (CLC    .LE.ZERO) CALL DDMEQNSUB (LSEQ(iK13),LSEN(jLC    ))
CCC      IF ( LMDRH ) THEN
CCC        IF (SC.EQ.'B' .OR. SC.EQ.'E' .OR. SC.EQ.'I') THEN
        IF (SC.EQ.'B' .OR. SC.EQ.'C' .OR.
     &      SC.EQ.'E' .OR. SC.EQ.'F' .OR.
     &      SC.EQ.'I' .OR. SC.EQ.'J') THEN ! ignore sensitivities of minor species
          CALL DDMEQNSUB (LSEQ(iMBCL ),LSEN(jCL  ))
          CALL DDMEQNSUB (LSEQ(iMBNO3),LSEN(jNO3 ))
          CALL DDMEQNSUB (LSEQ(iK3   ),LSEN(jHCL ))
          CALL DDMEQNSUB (LSEQ(iK4   ),LSEN(jHNO3))
          CALL DDMEQNSUB (LSEQ(iK2   ),LSEN(jNH3 ))
          LALFA = .TRUE.
        ENDIF
CCC      ENDIF
C     Get the dimension
      NDIM = 0
      DO I = 1, NSENS
        IF ( LSEQ(I) ) NDIM = NDIM + 1
      ENDDO
ccc      do ii = 1, nsens                            ! debug
ccc      if (lseq(ii)) write(*,*)'LSEQ:',ii,lseq(ii) ! debug
ccc      enddo                                       ! debug
C====================================================================
C     Fill the matrix A and constant vector B
C====================================================================
C
C     1. Equlibrium reactions
C
      CALL DDMDIFACT (DELGAMA)
      CX = 1.     ! Equlibrium equations are now normalized
C
C     K1
C
      iEQ = iK1
CCC   IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jH)/WATER*MOLAL(jSO4)/MOLAL(jHSO4) *
CX1     &       GAMA(mH2SO4)*(GAMA(mH2SO4)/GAMA(mHHSO4))**2.
        C1 =  3.*CX/GAMA(mH2SO4)
        C2 = -2.*CX/GAMA(mHHSO4)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mH2SO4) + C2*DELGAMA(K,mHHSO4)
        ENDDO
        AMAT(jH   ,iEQ) = AMAT(jH   ,iEQ) + CX/MAX(MOLAL(jH),TINY)
        AMAT(jSO4 ,iEQ) = AMAT(jSO4 ,iEQ) + CX/MAX(MOLAL(jSO4),TINY)
        AMAT(jHSO4,iEQ) = AMAT(jHSO4,iEQ) - CX/MAX(MOLAL(jHSO4),TINY)
        AMAT(jH2O ,iEQ) = AMAT(jH2O ,iEQ) - CX/WATER
CCC   ENDIF
C
C     K2
C
      iEQ = iK2
      IF ( LSEQ(iEQ) ) THEN
        IF (SC.EQ.'A') THEN
CX1        CX = MOLAL(jNH4)/MOLAL(jH)/GNH3 *
CX1     &       (GAMA(mNH4HS4)/GAMA(mHHSO4))**2. * XKW/R/TEMP
        C1 =  2.*CX/GAMA(mNH4HS4)
        C2 = -2.*CX/GAMA(mHHSO4)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNH4HS4) + C2*DELGAMA(K,mHHSO4)
        ENDDO
        ELSE
CX1        CX = MOLAL(jNH4)/MOLAL(jH)/GNH3 *
CX1     &       (GAMA(mNH4NO3)/GAMA(mHNO3))**2. * XKW/R/TEMP
        C1 =  2.*CX/GAMA(mNH4NO3)
        C2 = -2.*CX/GAMA(mHNO3)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNH4NO3) + C2*DELGAMA(K,mHNO3)
        ENDDO
        ENDIF
        AMAT(jNH4,iEQ) = AMAT(jNH4,iEQ) + CX/MAX(MOLAL(jNH4),TINY)
        AMAT(jH  ,iEQ) = AMAT(jH  ,iEQ) - CX/MAX(MOLAL(jH),TINY)
        AMAT(jNH3,iEQ) =                - CX/GNH3
      ENDIF
C
C     K3
C
      iEQ = iK3
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jH)/WATER*MOLAL(jCL)/WATER/GHCL *
CX1     &       GAMA(mHCL)**2. /R/TEMP
        C1 = 2.*CX/GAMA(mHCL)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mHCL)
        ENDDO
        AMAT(jH  ,iEQ) = AMAT(jH  ,iEQ) +    CX/MAX(MOLAL(jH),TINY)
        AMAT(jCL ,iEQ) = AMAT(jCL ,iEQ) +    CX/MAX(MOLAL(jCL),TINY)
        AMAT(jHCL,iEQ) =                -    CX/GHCL
        AMAT(jH2O,iEQ) = AMAT(jH2O,iEQ) - 2.*CX/WATER
      ENDIF
C
C     K4
C
      iEQ = iK4
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jH)/WATER*MOLAL(jNO3)/WATER/GHNO3 *
CX1     &       GAMA(mHNO3)**2. /R/TEMP
        C1 = 2.*CX/GAMA(mHNO3)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mHNO3)
        ENDDO
        AMAT(jH   ,iEQ) = AMAT(jH   ,iEQ) +    CX/MAX(MOLAL(jH),TINY)
        AMAT(jNO3 ,iEQ) = AMAT(jNO3 ,iEQ) +    CX/MAX(MOLAL(jNO3),TINY)
        AMAT(jHNO3,iEQ) =                 -    CX/GHNO3
        AMAT(jH2O ,iEQ) = AMAT(jH2O ,iEQ) - 2.*CX/WATER
      ENDIF
C
C     K5
C
      iEQ = iK5
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jNA)/WATER*MOLAL(jNA)/WATER*MOLAL(jSO4)/WATER *
CX1     &       GAMA(mNA2SO4)**3.
        C1 = 3.*CX/GAMA(mNA2SO4)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNA2SO4)
        ENDDO
        AMAT(jNA ,iEQ) = AMAT(jNA ,iEQ) + 2.*CX/MAX(MOLAL(jNA),TINY)
        AMAT(jSO4,iEQ) = AMAT(jSO4,iEQ) +    CX/MAX(MOLAL(jSO4),TINY)
        AMAT(jH2O,iEQ) = AMAT(jH2O,iEQ) - 3.*CX/WATER
      ENDIF
C
C     K6
C
      iEQ = iK6
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = GNH3*GHCL * (R*TEMP)**2.
        AMAT(jNH3,iEQ) = CX/GNH3
        AMAT(jHCL,iEQ) = CX/GHCL
      ENDIF
C
C     K7
C
      iEQ = iK7
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jNH4)/WATER*MOLAL(jNH4)/WATER*MOLAL(jSO4)/WATER *
CX1     &       GAMA(mNH42S4)**3.
        C1 = 3.*CX/GAMA(mNH42S4)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNH42S4)
        ENDDO
        AMAT(jNH4,iEQ) = AMAT(jNH4,iEQ) + 2.*CX/MAX(MOLAL(jNH4),TINY)
        AMAT(jSO4,iEQ) = AMAT(jSO4,iEQ) +    CX/MAX(MOLAL(jSO4),TINY)
        AMAT(jH2O,iEQ) = AMAT(jH2O,iEQ) - 3.*CX/WATER
      ENDIF
C
C     K8
C
      iEQ = iK8
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jNA)/WATER*MOLAL(jCL)/WATER *
CX1     &       GAMA(mNACL)**2.
        C1 = 2.*CX/GAMA(mNACL)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNACL)
        ENDDO
        AMAT(jNA ,iEQ) = AMAT(jNA ,iEQ) +    CX/MAX(MOLAL(jNA),TINY)
        AMAT(jCL ,iEQ) = AMAT(jCL ,iEQ) +    CX/MAX(MOLAL(jCL),TINY)
        AMAT(jH2O,iEQ) = AMAT(jH2O,iEQ) - 2.*CX/WATER
      ENDIF
C
C     K9
C
      iEQ = iK9
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jNA)/WATER*MOLAL(jNO3)/WATER *
CX1     &       GAMA(mNANO3)**2.
        C1 = 2.*CX/GAMA(mNANO3)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNANO3)
        ENDDO
        AMAT(jNA ,iEQ) = AMAT(jNA ,iEQ) +    CX/MAX(MOLAL(jNA),TINY)
        AMAT(jNO3,iEQ) = AMAT(jNO3,iEQ) +    CX/MAX(MOLAL(jNO3),TINY)
        AMAT(jH2O,iEQ) = AMAT(jH2O,iEQ) - 2.*CX/WATER
      ENDIF
C
C     K10
C
      iEQ = iK10
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = GNH3*GHNO3 * (R*TEMP)**2.
        AMAT(jNH3 ,iEQ) = CX/GNH3
        AMAT(jHNO3,iEQ) = CX/GHNO3
      ENDIF
C
C     K11
C
      iEQ = iK11
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jNA)/WATER*MOLAL(jHSO4)/WATER *
CX1     &       GAMA(mNAHSO4)**2.
        C1 = 2.*CX/GAMA(mNAHSO4)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNAHSO4)
        ENDDO
        AMAT(jNA  ,iEQ) = AMAT(jNA  ,iEQ) +    CX/MAX(MOLAL(jNA),TINY)
        AMAT(jHSO4,iEQ) = AMAT(jHSO4,iEQ) +    CX/MAX(MOLAL(jHSO4),TINY)
        AMAT(jH2O ,iEQ) = AMAT(jH2O ,iEQ) - 2.*CX/WATER
      ENDIF
C
C     K12
C
      iEQ = iK12
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = MOLAL(jNH4)/WATER*MOLAL(jHSO4)/WATER *
CX1     &       GAMA(mNH4HS4)**2.
        C1 = 2.*CX/GAMA(mNH4HS4)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mNH4HS4)
        ENDDO
        AMAT(jNH4 ,iEQ) = AMAT(jNH4 ,iEQ) +    CX/MAX(MOLAL(jNH4),TINY)
        AMAT(jHSO4,iEQ) = AMAT(jHSO4,iEQ) +    CX/MAX(MOLAL(jHSO4),TINY)
        AMAT(jH2O ,iEQ) = AMAT(jH2O ,iEQ) - 2.*CX/WATER
      ENDIF
C
C     K13
C
      iEQ = iK13
      IF ( LSEQ(iEQ) ) THEN
CX1        CX = (MOLAL(jNH4)/WATER)**3. * 
CX1     &       (MOLAL(jHSO4)/WATER*MOLAL(jSO4)/WATER) *
CX1     &       GAMA(mLC)**5.
        C1 = 5.*CX/GAMA(mLC)
        DO K = 1, NIONSP1
          AMAT(K,iEQ) = C1*DELGAMA(K,mLC)
        ENDDO
        AMAT(jNH4 ,iEQ) = AMAT(jNH4 ,iEQ) + 3.*CX/MAX(MOLAL(jNH4),TINY)
        AMAT(jHSO4,iEQ) = AMAT(jHSO4,iEQ) +    CX/MAX(MOLAL(jHSO4),TINY)
        AMAT(jSO4 ,iEQ) = AMAT(jSO4 ,iEQ) +    CX/MAX(MOLAL(jSO4),TINY)
        AMAT(jH2O ,iEQ) = AMAT(jH2O ,iEQ) - 5.*CX/WATER
      ENDIF
C
C     2. Mass balances
C
      iEQ = iMBNA
      IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XLIQ(      iEQ,J) = SINI(jTNA,J)
        ENDDO
        AMAT(jNA    ,iEQ) = 1.
        AMAT(jNACL  ,iEQ) = 1.
        AMAT(jNA2SO4,iEQ) = 2.
        AMAT(jNANO3 ,iEQ) = 1.
        AMAT(jNAHSO4,iEQ) = 1.
      ENDIF

      iEQ = iMBSO4
CCC   IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XLIQ(      iEQ,J) = SINI(jTSO4,J)
        ENDDO
        AMAT(jSO4   ,iEQ) = 1.
        AMAT(jHSO4  ,iEQ) = 1.
        AMAT(jNH42S4,iEQ) = 1.
        AMAT(jNH4HS4,iEQ) = 1.
        AMAT(jNA2SO4,iEQ) = 1.
        AMAT(jNAHSO4,iEQ) = 1.
        AMAT(jLC    ,iEQ) = 2.
CCC   ENDIF

      iEQ = iMBNH4
CCC   IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XLIQ(      iEQ,J) = SINI(jTNH4,J)
        ENDDO
        AMAT(jNH3   ,iEQ) = 1.
        AMAT(jNH4   ,iEQ) = 1.
        AMAT(jNH42S4,iEQ) = 2.
        AMAT(jNH4HS4,iEQ) = 1.
        AMAT(jNH4NO3,iEQ) = 1.
        AMAT(jNH4CL ,iEQ) = 1.
        AMAT(jLC    ,iEQ) = 3.
CCC   ENDIF

      iEQ = iMBNO3
      IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XLIQ(      iEQ,J) = SINI(jTNO3,J)
        ENDDO
        AMAT(jHNO3  ,iEQ) = 1.
        AMAT(jNO3   ,iEQ) = 1.
        AMAT(jNANO3 ,iEQ) = 1.
        AMAT(jNH4NO3,iEQ) = 1.
      ENDIF

      iEQ = iMBCL
      IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          XLIQ(      iEQ,J) = SINI(jTCL,J)
        ENDDO
        AMAT(jHCL   ,iEQ) = 1.
        AMAT(jCL    ,iEQ) = 1.
        AMAT(jNACL  ,iEQ) = 1.
        AMAT(jNH4CL ,iEQ) = 1.
      ENDIF
C
C     3. Charge balances
C
      iEQ = iCB
CCC   IF ( LSEQ(iEQ) ) THEN
        AMAT(jH   ,iEQ) =  1. + XKW*RH*(WATER/MAX(MOLAL(jH),TINY))**2.
        AMAT(jNA  ,iEQ) =  1.
        AMAT(jNH4 ,iEQ) =  1.
        AMAT(jSO4 ,iEQ) = -2.
        AMAT(jHSO4,iEQ) = -1.
        AMAT(jNO3 ,iEQ) = -1.
        AMAT(jCL  ,iEQ) = -1.
        AMAT(jH2O ,iEQ) = -2.*XKW*RH*WATER/MAX(MOLAL(jH),TINY)
CCC   ENDIF
C
C     4. ZSR relationship
C
      iEQ = iZSR
CCC   IF ( LSEQ(iEQ) ) THEN
        DO J = 1, NDDMSP
          CALL DDMDIFZSR (SINI(1,J),AMAT(1,iEQ),XLIQ(iEQ,J),IERR,IOUT)
          IF (IERR.LT.0) RETURN
        ENDDO
CCC   ENDIF
C====================================================================
C
C     Solve the sensitivity eqns
C
      CALL DDMEQNSLV (NSENS,LSEQ,LSEN,AMAT,NDDMSP,XLIQ,
     &                NDIM,AA,BB,IPVT,IERR,IOUT)
      IF (IERR.LT.0) THEN
        WRITE(IOUT,'(A)') ' Called by DDMISRLIQ'
        RETURN
      ENDIF
C
C     Calculate DALFA for minor species equilibria
C
      IF ( LALFA ) THEN
CCC     A1 = XK3*R*TEMP*(WATER/GAMA(mHCL))**2.0
CCC     A2 = XK4*R*TEMP*(WATER/GAMA(mHNO3))**2.0
CCC     A3 = (XKW/XK2)/(R*TEMP)*(GAMA(mNH4NO3)/GAMA(mHNO3))**2.0
CCC     DO J = 1, NDDMSP
CCC       XSUM1 = 0.0
CCC       XSUM2 = 0.0
CCC       XSUM3 = 0.0
CCC       DO K = 1, NIONSP1
CCC         XSUM1 = XSUM1 + DELGAMA(K,mHCL   )*XLIQ(K,J)
CCC         XSUM2 = XSUM2 + DELGAMA(K,mHNO3  )*XLIQ(K,J)
CCC         XSUM3 = XSUM3 + DELGAMA(K,mNH4NO3)*XLIQ(K,J)
CCC       ENDDO
CCC       DALFA(1,J) = 2.*A1*( XLIQ(jH2O,J)/WATER -XSUM1/GAMA(mHCL ) )
CCC       DALFA(2,J) = 2.*A2*( XLIQ(jH2O,J)/WATER -XSUM2/GAMA(mHNO3) )
CCC       DALFA(3,J) = 2.*A3*( XSUM3/GAMA(mNH4NO3)-XSUM2/GAMA(mHNO3) )
CCC     ENDDO
CCC   This is not always working - Let's make it simple
        DO J = 1, NDDMSP
          XLIQ(jHCL ,J) = SINI(jTCL ,J)
          XLIQ(jHNO3,J) = SINI(jTNO3,J)
          XLIQ(jNH3 ,J) = 0.0
        ENDDO
      ENDIF

      RETURN
      END
C     END OF SUBROUTINE DDMISRLIQ


      SUBROUTINE DDMDIFACT (DELGAMA)

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL*8    LN10                      ! Natural log of 10 (see CALCACT)
      PARAMETER (LN10=2.30258509299404568402D0)

      REAL      DELI(NIONSP1)             ! del IONIC by del Y(A)
      REAL      G0(NPAIR)                 ! Common log of binary activity coefficients
      REAL      DELG0(NIONSP1,NPAIR)      ! del G0 by del Y(A)

      REAL      ZPL,ZMI,AGAMA,SION,H,CH,HD
      REAL*8    XPL,XMI,XIJ,YJI,DELXIJ,DELYJI
      REAL      DELF1(NIONSP1,3),DELF2(NIONSP1,4)

      REAL      DELGAMA(NIONSP1,NPAIR)    ! del GAMA by del Y(A)
      INTEGER   IJMAP(3,4),I,J,K
C
C     Ionic strength
C
      IF (IONIC.GE.100.d0) THEN           ! IONIC is limited in (TINY, 100.d0) in ISORROPIA
        DO I = 1, NIONSP1
          DELI(I) = 0.0
        ENDDO
      ELSE
        DELI(jH2O) = 0.0
        DO I = 1, NIONS
          DELI(I) = 0.5*Z(I)*Z(I)/WATER
          DELI(jH2O) = DELI(jH2O) + MOLAL(I)*Z(I)*Z(I)
        ENDDO
        DELI(jH2O) = -0.5*DELI(jH2O)/(WATER*WATER)
      ENDIF
C
C     Binary activity coefficients (Kusik and Meissner method)
C
      CALL DDMKMFUL (NIONSP1,NPAIR,IONIC,SNGL(TEMP),DELI,G0,DELG0)
C
C     Mapping of IJMAP
C
      IJMAP(1,1) = mHCL
      IJMAP(1,2) = mH2SO4
      IJMAP(1,3) = mHHSO4
      IJMAP(1,4) = mHNO3
      IJMAP(2,1) = mNACL
      IJMAP(2,2) = mNA2SO4
      IJMAP(2,3) = mNAHSO4
      IJMAP(2,4) = mNANO3
      IJMAP(3,1) = mNH4CL
      IJMAP(3,2) = mNH42S4
      IJMAP(3,3) = mNH4HS4
      IJMAP(3,4) = mNH4NO3
C
C     Multicomponent activity coefficients (Bromley's method)
C
      AGAMA = 0.511*(298.0/TEMP)**1.5     ! Debye Huckel const. at T
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1+SION)
      HD    = 0.5*AGAMA/(SION*(1.+SION)*(1.+SION))

      DO K = 1, NIONSP1
        DO 100 I=1,3
           DELF1(K,I)=0.0
           DELF2(K,I)=0.0
100     CONTINUE
        DELF2(K,4)=0.0
      ENDDO

      DO 110 I=1,3
         ZPL = Z(I)
         XPL = MOLAL(I)/WATER
         DO 110 J=1,4
            ZMI   = Z(J+3)
            XMI   = MOLAL(J+3)/WATER
            CH    = 0.25*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*XPL
            YJI   = CH*XMI

            DO K = 1, NIONSP1
              DELXIJ = -XPL*DELI(K)/IONIC
              DELYJI = -XMI*DELI(K)/IONIC
              IF (K.EQ.I) THEN
                DELXIJ = DELXIJ + 1./WATER
              ELSEIF (K.EQ.J+3) THEN
                DELYJI = DELYJI + 1./WATER
              ELSEIF (K.EQ.jH2O) THEN
                DELXIJ = DELXIJ - XPL/WATER
                DELYJI = DELYJI - XMI/WATER
              ENDIF
              DELXIJ = CH*DELXIJ
              DELYJI = CH*DELYJI

              DELF1(K,I) = DELF1(K,I)
     &                   + G0(IJMAP(I,J)) * DELYJI
     &                   + YJI * DELG0(K,IJMAP(I,J))
     &                   + ZPL*ZMI*HD * YJI * DELI(K)
     &                   + ZPL*ZMI*H * DELYJI
              DELF2(K,J) = DELF2(K,J)
     &                   + G0(IJMAP(I,J)) * DELXIJ
     &                   + XIJ * DELG0(K,IJMAP(I,J))
     &                   + ZPL*ZMI*HD * XIJ * DELI(K)
     &                   + ZPL*ZMI*H * DELXIJ
            ENDDO

110   CONTINUE
C
C     del log10(GAMA)
C
      DO 120 I=1,3
         ZPL = Z(I)
         DO 120 J=1,4
            ZMI = Z(J+3)

            DO K = 1, NIONSP1
               DELGAMA(K,IJMAP(I,J)) = ZPL*ZMI * (
     &                  (DELF1(K,I)/ZPL + DELF2(K,J)/ZMI) / (ZPL+ZMI)
     &                               - HD * DELI(K) )
            ENDDO

120   CONTINUE

      DO K = 1, NIONSP1
         DELGAMA(K,mLC) = 0.20 * ( 3.0*DELGAMA(K,mNH42S4)
     &                           + 2.0*DELGAMA(K,mNH4HS4) )
      ENDDO
C
C     Convert del log10(GAMA) to del GAMA
C
      DO I = 1, NPAIR
         IF (GAMA(I).LE.1.d-5 .OR. GAMA(I).GE.1.d5) THEN
            DO K = 1, NIONSP1
               DELGAMA(K,I) = 0.0
            ENDDO
         ELSE
            DO K = 1, NIONSP1
               DELGAMA(K,I) = LN10 * GAMA(I) * DELGAMA(K,I)
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMDIFACT


      SUBROUTINE DDMKMFUL (N,NPAIRS,IONIC,TEMP,DELI,G0,DELG0)

      INTEGER   N, NPAIRS
      REAL      IONIC, TEMP, DELI(N), G0(NPAIRS), DELG0(N,NPAIRS)
      REAL      SION, CUBI, TI, CF1, CF2, CF2D

      INTEGER   NPAIRD, I, K
      PARAMETER (NPAIRD=10)               ! Number of ion pairs whose Q value is available
      INTEGER   IG(NPAIRD)
      DATA IG / 1,2,3,4,5,6,7,8,10,11 /
      REAL      ZI(NPAIRD)                ! Mapping of Q to the internal order of ion pairs
      DATA ZI / 1., 2., 1., 2., 1., 1., 2., 1., 1., 1. /
      REAL      Q(NPAIRD)                 ! Kusik-Meissner parameters (see KMFUL)
      DATA Q  / 2.23,-0.19,-0.39,-0.25,-1.15,0.82,-0.1,8.0,2.6,6.0 /
C
      SION = SQRT(IONIC)
      CUBI = IONIC*IONIC*IONIC
C
C     Coefficients at 25 oC
C
      DO I = 1, NPAIRD
         CALL DDMMKBI(N,Q(I),IONIC,SION,CUBI,ZI(I),G0(IG(I)),DELI,
     &                                          DELG0(1,IG(I)))
      ENDDO
C
C     Correct for T other than 298 K
C
      TI  = TEMP-273.0
      IF (ABS(TI-25.0) .GT. 1.0) THEN
         CF1 = 1.125-0.005*TI
         CF2 = (CF1-1.)*(0.039*IONIC**0.92-0.41*SION/(1.+SION))
         CF2D = (CF1-1.)*( .03588/IONIC**.08
     &                    -.205/(SION*(1.+SION)*(1.+SION)) )
         DO I = 1, NPAIRD
            G0(IG(I)) = CF1*G0(IG(I)) - CF2*ZI(I)
            DO K = 1, N
               DELG0(K,IG(I)) = CF1*DELG0(K,IG(I)) - ZI(I)*CF2D*DELI(K)
            ENDDO
         ENDDO
      ENDIF
C
      G0( 9) = G0( 6) + G0( 8) - G0(11)
      G0(12) = G0( 1) + G0( 8) - G0(11)
      DO K = 1, N
         DELG0(K, 9) = DELG0(K, 6) + DELG0(K, 8) - DELG0(K,11)
         DELG0(K,12) = DELG0(K, 1) + DELG0(K, 8) - DELG0(K,11)
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMKMFUL


      SUBROUTINE DDMMKBI (N,Q,IONIC,SION,CUBI,ZIP,G,DELI,DELG)
C
      REAL      LN10                      ! Natural log of 10 (see INIT1,INIT2,INIT3)
      PARAMETER (LN10=2.3025851)
      INTEGER   N
      REAL      Q, IONIC, SION, CUBI, ZIP, G, DELI(N), DELG(N)

      REAL      B, C, XX, BI, XX1, XX2, XX3, XX4
      INTEGER   K
C
      B=.75-.065*Q
      C= 1.0
      IF (IONIC.LT.6.0) C=1.+.055*Q*EXP(-.023*CUBI)
      XX=-0.5107*SION/(1.+C*SION)
      BI=(1.+B*(1.+.1*IONIC)**Q-B)
      G =ZIP*ALOG10(BI) + ZIP*XX

      XX1 = .1*B*Q*(1.+.1*IONIC)**(Q-1.)/(BI*LN10)
      XX2 = 0.5/SION+.003795*Q*CUBI*EXP(-.023*CUBI)
      XX3 = (1.+C*SION)*(1.+C*SION)
      XX4 = ZIP*(XX1-.5107*XX2/XX3)
      DO K = 1, N
         DELG(K) = XX4 * DELI(K)
      ENDDO
C
      RETURN
      END
C     END OF SUBROUTINE DDMMKBI


      SUBROUTINE DDMDIFZSR (SI,AVEC,CONST,IERR,IOUT)

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL      SI(NINPT)                 ! Initial Sensitivities of total
      REAL      AVEC(NSENS), CONST
      REAL      TVEC(NSENS)               ! Temporary solution vector for case I
      INTEGER   IERR, IOUT

      CHARACTER SC*1
      REAL*8    SO4I, HSO4I, AML5, FRNH4, FRNA, FRNO3, FRCL
      REAL      DFRNA
      INTEGER   I

      AVEC(jH2O) = 1.

      SC = SCASE(1:1)
      IF     (SC.EQ.'A') THEN
        AVEC(jSO4   ) = -1./M0(mNH42S4)
        AVEC(jHSO4  ) = -1./M0(mNH42S4)
      ELSEIF (SC.EQ.'C' .OR. SC.EQ.'F' .OR. SC.EQ.'J') THEN
        AVEC(jSO4   ) = -1./M0(mH2SO4 )
        AVEC(jHSO4  ) = -1./M0(mH2SO4 )
        AVEC(jNH4   ) =  1./M0(mH2SO4 ) -1./M0(mNH4HS4)
        AVEC(jNA    ) =  1./M0(mH2SO4 ) -1./M0(mNAHSO4)
      ELSEIF (SC.EQ.'B' .OR. SC.EQ.'E') THEN
        SO4I  = MOLAL(jSO4 ) -MOLAL(jH)
        HSO4I = MOLAL(jHSO4) +MOLAL(jH)
        IF (SO4I.LT.HSO4I) THEN
          AVEC(jHSO4) = -1./M0(mNH4HS4)
          AVEC(jSO4 ) =  1./M0(mNH4HS4) -1./M0(mLC)
          AVEC(jH   ) = -2./M0(mNH4HS4) +1./M0(mLC)
        ELSE
          AVEC(jSO4 ) = -1./M0(mNH42S4)
          AVEC(jHSO4) =  1./M0(mNH42S4) -1./M0(mLC)
          AVEC(jH   ) =  2./M0(mNH42S4) -1./M0(mLC)
        ENDIF
      ELSEIF (SC.EQ.'D') THEN
        AVEC(jSO4   ) = -1./M0(mNH42S4)
        AVEC(jHSO4  ) = -1./M0(mNH42S4)
        AML5 = MOLAL(jNH4)-2.D0*(MOLAL(jSO4)+MOLAL(jHSO4))
        IF (AML5.LT.MOLAL(jNO3)) THEN
          AVEC(jNH4 ) =             -1./M0(mNH4NO3)
          AVEC(jSO4 ) = AVEC(jSO4 ) +2./M0(mNH4NO3)
          AVEC(jHSO4) = AVEC(jHSO4) +2./M0(mNH4NO3)
        ELSE
          AVEC(jNO3 ) =             -1./M0(mNH4NO3)
        ENDIF
      ELSEIF (SC.EQ.'G') THEN
        AVEC(jSO4   ) = -1./M0(mNH42S4)
        AVEC(jHSO4  ) = -1./M0(mNH42S4)
        AVEC(jNA    ) =  .5/M0(mNH42S4) -.5/M0(mNA2SO4)
        FRNH4 = MOLAL(jNH4)+MOLAL(jNA)-2.D0*(MOLAL(jSO4)+MOLAL(jHSO4))
        IF (FRNH4.LT.MOLAL(jNO3)) THEN
          AVEC(jNH4 ) =               -1./M0(mNH4NO3)
          AVEC(jNA  ) = AVEC(jNA  )   -1./M0(mNH4NO3)
          AVEC(jSO4 ) = AVEC(jSO4 )   +2./M0(mNH4NO3)
          AVEC(jHSO4) = AVEC(jHSO4)   +2./M0(mNH4NO3)
        ELSE
          AVEC(jNO3 ) =               -1./M0(mNH4NO3)
          FRNH4 = FRNH4 -MOLAL(jNO3)
          IF (FRNH4.LT.MOLAL(jCL)) THEN
            AVEC(jNH4 ) =             -1./M0(mNH4CL)
            AVEC(jNA  ) = AVEC(jNA  ) -1./M0(mNH4CL)
            AVEC(jSO4 ) = AVEC(jSO4 ) +2./M0(mNH4CL)
            AVEC(jHSO4) = AVEC(jHSO4) +2./M0(mNH4CL)
            AVEC(jNO3 ) = AVEC(jNO3 ) +1./M0(mNH4CL)
          ELSE
            AVEC(jCL  ) =             -1./M0(mNH4CL)
          ENDIF
        ENDIF
      ELSEIF (SC.EQ.'H') THEN
        CONST   = SI(jTSO4)/M0(mNA2SO4)
        AVEC(jNA2SO4) =  1./M0(mNA2SO4)
        FRNA = W(jTNA) -2.D0*W(jTSO4)
        DFRNA=SI(jTNA) -2. *SI(jTSO4)
        IF (FRNA.LT.W(jTNO3)) THEN
          CONST = CONST +DFRNA/M0(mNANO3)
          AVEC(jNANO3) =    1./M0(mNANO3)
          FRNO3 = MOLAL(jNO3) -FRNA +CNANO3
          IF (FRNO3.GT.1.D2*TINY) THEN
            IF (FRNO3.LT.MOLAL(jNH4)) THEN
              CONST = CONST              -DFRNA/M0(mNH4NO3)
              AVEC(jNO3  ) =                -1./M0(mNH4NO3)
              AVEC(jNANO3) = AVEC(jNANO3)   -1./M0(mNH4NO3)
              FRNH4 = MOLAL(jNH4) -FRNO3
              IF (FRNH4.LT.MOLAL(jCL)) THEN
                CONST = CONST            +DFRNA/M0(mNH4CL)
                AVEC(jNH4  ) =              -1./M0(mNH4CL)
                AVEC(jNO3  ) = AVEC(jNO3  ) +1./M0(mNH4CL)
                AVEC(jNANO3) = AVEC(jNANO3) +1./M0(mNH4CL)
              ELSE
                AVEC(jCL   ) =              -1./M0(mNH4CL)
              ENDIF
            ELSE
              AVEC(jNH4  ) =                -1./M0(mNH4NO3)
            ENDIF
          ELSE
            FRNH4 = MOLAL(jNH4) -FRNO3
            IF (FRNH4.GT.ZERO) THEN
              IF (MOLAL(jCL).LT.FRNH4) THEN
                AVEC(jCL   ) =              -1./M0(mNH4CL)
              ELSE
                AVEC(jNH4  ) =              -1./M0(mNH4CL)
              ENDIF
            ENDIF
          ENDIF
        ELSE
          CONST = CONST +SI(jTNO3)/M0(mNANO3)
          AVEC(jNANO3) =        1./M0(mNANO3)
          FRNA = FRNA -W(jTNO3)
          DFRNA=DFRNA-SI(jTNO3)
          IF (FRNA.LT.W(jTCL)) THEN
            CONST = CONST +DFRNA/M0(mNACL)
            AVEC(jNACL) =     1./M0(mNACL)
            FRCL = MOLAL(jCL) -FRNA +CNACL
            FRNH4= MOLAL(jNH4) -TINY
            IF (FRNH4.GT.ZERO) THEN
              IF (FRCL.LT.FRNH4) THEN
                CONST = CONST            -DFRNA/M0(mNH4CL)
                AVEC(jCL  ) =               -1./M0(mNH4CL)
                AVEC(jNACL) = AVEC(jNACL)   -1./M0(mNH4CL)
              ELSE
                AVEC(jNH4  ) =              -1./M0(mNH4CL)
              ENDIF
            ENDIF
          ELSE
            CONST = CONST +SI(jTCL)/M0(mNACL)
            AVEC(jNACL) =        1./M0(mNACL)
          ENDIF
        ENDIF
      ELSEIF (SC.EQ.'I') THEN
        ! Get dY/dp from solid-only case
        DO I = 1, NSENS
          TVEC(I) = 0.0
        ENDDO
        CALL DDMISRDRY (SI,TVEC,IERR,IOUT)
        IF (IERR.LT.0) RETURN
        ! dE/dp = dY/dp(from solid-only case) - dY/dp(from current case)
        CONST = TVEC(jNA2SO4)/M0(mNA2SO4) +TVEC(jNAHSO4)/M0(mNAHSO4)
     &         +TVEC(jNH42S4)/M0(mNH42S4) +TVEC(jNH4HS4)/M0(mNH4HS4)
     &         +TVEC(jLC)    /M0(mLC)
        AVEC(jNA2SO4) = 1./M0(mNA2SO4)
        AVEC(jNAHSO4) = 1./M0(mNAHSO4)
        AVEC(jNH42S4) = 1./M0(mNH42S4)
        AVEC(jNH4HS4) = 1./M0(mNH4HS4)
        AVEC(jLC)     = 1./M0(mLC)
      ELSE
        WRITE(IOUT,'(//,A)') 'ERROR in DDMDIFZSR:'
        WRITE(IOUT,'(A,A)')  ' Case not supported - ',SCASE
        IERR = -1
        RETURN
      ENDIF

      RETURN
      END
C     END OF SUBROUTINE DDMDIFZSR


      SUBROUTINE DDMISRMDR (CMDRH,WF)

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL*8    WF
      CHARACTER CMDRH*1                   ! 'D' - dry case; 'L' - saturated liquid case

      REAL      XSAT(NSENS,MXFDDM)        ! Saturated liquid solution of MDRH region
      REAL      DDAMSUL, DDSOSUL, DDAMBIS, DDSOBIS, DDLC, DDAMNIT,
     &          DDAMCHL, DDSONIT, DDSOCHL, DDAMG, DDHAG, DDNAG
      REAL*8    ONEMWF
      INTEGER   IERR, IOUT, I, J
C
C     Initial check
C
      IF ( .NOT.LDDM ) RETURN             ! Check if DDM is being used
      IF ( IPAR.LT.0 ) RETURN             ! Check error from the previous call
      IF ( SCASE.EQ.'G2 ; SUBCASE 3' ) RETURN   ! Skip [G2 ; SUBCASE 3]
C
C     Initialize error flag & unit number of log file
C
      IERR = 1
      IOUT = IPAR
C
C     Solutions of 'dry' and 'saturated liquid' cases
C
      IF (CMDRH.EQ.'D') THEN              ! Dry case
        LMDRH = .TRUE.                    ! Set LMDRH flag
        DO J = 1, NDDMSP
          CALL DDMISRDRY (SINI(1,J),XMAT(1,J),IERR,IOUT)
          IF (IERR.LT.0) GOTO 900
        ENDDO
      ELSE                                ! Saturated liquid case
        IF ( IONIC.LT.TINY2 ) RETURN
        DO J = 1, NDDMSP
          DO I = 1, NSENS
            XSAT(I,J) = 0.0               ! Initialize the solution matrix
          ENDDO
        ENDDO
CCC        MOLAL(jH     ) = MAX(MOLAL(jH     ), TINY)
CCC        MOLAL(jHSO4  ) = MAX(MOLAL(jHSO4  ), TINY)
        CALL DDMISRLIQ (XSAT,IERR,IOUT)
        IF (IERR.LT.0) GOTO 900
C
C     Sum of two weighted solutions
C
        ONEMWF  = ONE - WF
        IF (SCASE(1:1).EQ.'B') THEN ! cbk=> ( MOLAL(5) and MOLAL(6) are differently assigned 
          DO J = 1, NDDMSP
            XSAT(jHSO4,J) = 0.0
          ENDDO
        ENDIF                       ! cbk<=   in CALCB2A2 and CALCB1B than in CALCMDRH. Why? )
        DO J = 1, NDDMSP
C     Sensitivities of salt dissolutions
          DDAMSUL = XMAT(jNH42S4,J) - XSAT(jNH42S4,J)
          DDSOSUL = XMAT(jNA2SO4,J) - XSAT(jNA2SO4,J)
          DDAMBIS = XMAT(jNH4HS4,J) - XSAT(jNH4HS4,J)
          DDSOBIS = XMAT(jNAHSO4,J) - XSAT(jNAHSO4,J)
          DDLC    = XMAT(jLC    ,J) - XSAT(jLC    ,J)
          DDAMNIT = XMAT(jNH4NO3,J) - XSAT(jNH4NO3,J)
          DDAMCHL = XMAT(jNH4CL ,J) - XSAT(jNH4CL ,J)
          DDSONIT = XMAT(jNANO3 ,J) - XSAT(jNANO3 ,J)
          DDSOCHL = XMAT(jNACL  ,J) - XSAT(jNACL  ,J)
C     Sensitivities of gas dissolutions
          DDAMG   = XMAT(jNH3   ,J) - XSAT(jNH3   ,J)
          DDHAG   = XMAT(jHCL   ,J) - XSAT(jHCL   ,J)
          DDNAG   = XMAT(jHNO3  ,J) - XSAT(jHNO3  ,J)
C     Sensitivities at MDRH
          XMAT(jH   ,J) = ONEMWF*XSAT(jH   ,J)
          XMAT(jNA  ,J) = ONEMWF*( 2.*DDSOSUL +DDSOBIS +DDSONIT
     &                                                 +DDSOCHL )
          XMAT(jNH4 ,J) = ONEMWF*( 2.*DDAMSUL +DDAMBIS +3.*DDLC
     &                               +DDAMNIT +DDAMCHL +DDAMG   )
          XMAT(jCL  ,J) = ONEMWF*(    DDAMCHL +DDSOCHL +DDHAG   )
          XMAT(jSO4 ,J) = ONEMWF*(    DDAMSUL +DDSOSUL +DDLC
     &                               -XSAT(jHSO4,J)             )
          XMAT(jHSO4,J) = ONEMWF*(    DDAMBIS +DDSOBIS +DDLC
     &                               +XSAT(jHSO4,J)             )
          XMAT(jNO3 ,J) = ONEMWF*(    DDAMNIT +DDSONIT +DDNAG   )
          XMAT(jH2O ,J) = ONEMWF*XSAT(jH2O ,J)
          DO I = jH2O + 1, NSENS
            XMAT(I,J) = WF*XMAT(I,J) +ONEMWF*XSAT(I,J)
          ENDDO
        ENDDO
      ENDIF

      RETURN

900   CONTINUE
      WRITE(IOUT,'(A)')' The above error occurred in the MDRH region'
      IPAR = IERR
      RETURN

      END
C     END OF SUBROUTINE DDMISRMDR


      SUBROUTINE DDMCALCHA

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL*8    ALFA
      INTEGER   J

      IF (WATER.GT.TINY) THEN
        ALFA = XK3*R*TEMP*(WATER/GAMA(mHCL))**2.0
        DO J = 1, NDDMSP
          XMAT(jCL  ,J) = ( SINI(jTCL,J)*ALFA -XMAT(jH,J)*MOLAL(jCL)
     &                     +GHCL*DALFA(1,J) ) /
     &                    ( MOLAL(jCL) +MOLAL(jH) +ALFA )
        ENDDO
      ENDIF
      DO J = 1, NDDMSP
        XMAT(jHCL ,J) = SINI(jTCL ,J) -XMAT(jCL,J)
        XMAT(jH   ,J) = XMAT(jH   ,J) +XMAT(jCL,J)
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMCALCHA


      SUBROUTINE DDMCALCNA

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL*8    ALFA
      INTEGER   J

      IF (WATER.GT.TINY) THEN
        ALFA = XK4*R*TEMP*(WATER/GAMA(mHNO3))**2.0
        DO J = 1, NDDMSP
          XMAT(jNO3 ,J) = ( SINI(jTNO3,J)*ALFA -XMAT(jH,J)*MOLAL(jNO3)
     &                     +GHNO3*DALFA(2,J) ) /
     &                    ( MOLAL(jNO3) +MOLAL(jH) +ALFA )
        ENDDO
      ENDIF
      DO J = 1, NDDMSP
        XMAT(jHNO3,J) = SINI(jTNO3,J) -XMAT(jNO3,J)
        XMAT(jH   ,J) = XMAT(jH   ,J) +XMAT(jNO3,J)
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMCALCNA


      SUBROUTINE DDMCALCNHA

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL*8    ALFA1, ALFA2, X2NU, X2DE
      INTEGER   J

      IF (W(jTCL).LE.TINY .AND. W(jTNO3).LE.TINY) THEN
        RETURN
      ELSEIF (W(jTCL).LE.TINY) THEN
        CALL DDMCALCNA
        RETURN
      ELSEIF (W(jTNO3).LE.TINY) THEN
        CALL DDMCALCHA
        RETURN
      ENDIF

      IF (WATER.GT.TINY) THEN
        ALFA1 = XK3*R*TEMP*(WATER/GAMA(mHCL ))**2.0 ! HCL
        ALFA2 = XK4*R*TEMP*(WATER/GAMA(mHNO3))**2.0 ! HNO3
        DO J = 1, NDDMSP
          X2NU = SINI(jTNO3,J)*ALFA2 -XMAT(jH,J)*MOLAL(jNO3)
     &          +GHNO3*DALFA(2,J)
          X2DE = MOLAL(jNO3) +MOLAL(jH) +ALFA2
          XMAT(jCL  ,J) = ( SINI(jTCL,J)*ALFA1 +GHCL*DALFA(1,J)
     &                     -MOLAL(jCL)*(XMAT(jH,J) +X2NU/X2DE) ) /
     &                    ( MOLAL(jCL) +MOLAL(jH) +ALFA1
     &                     -MOLAL(jCL)*MOLAL(jNO3)/X2DE )
          XMAT(jNO3 ,J) = ( X2NU -MOLAL(jNO3)*XMAT(jCL,J) ) / X2DE
        ENDDO
      ENDIF
      DO J = 1, NDDMSP
        XMAT(jHCL ,J) = SINI(jTCL ,J) -XMAT(jCL ,J)
        XMAT(jHNO3,J) = SINI(jTNO3,J) -XMAT(jNO3,J)
        XMAT(jH   ,J) = XMAT(jH   ,J) +XMAT(jCL ,J) +XMAT(jNO3,J)
      ENDDO

      RETURN
      END
C     END OF SUBROUTINE DDMCALCNHA


      SUBROUTINE DDMCALCNH3

      use tracer
      INCLUDE 'isrpia.inc'
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'

      REAL*8    ALFA
      INTEGER   J

      IF (WATER.GT.TINY) THEN
        ALFA = (XKW/XK2)/(R*TEMP)*(GAMA(mNH4NO3)/GAMA(mHNO3))**2.0
        DO J = 1, NDDMSP
          XMAT(jNH3 ,J) = ( XMAT(jNH4,J)*ALFA -XMAT(jH,J)*GNH3
     &                     +MOLAL(jNH4)*DALFA(3,J) ) /
     &                    ( GNH3 +MOLAL(jH) +ALFA )
          XMAT(jNH4 ,J) = XMAT(jNH4 ,J) -XMAT(jNH3,J)
          XMAT(jH   ,J) = XMAT(jH   ,J) +XMAT(jNH3,J)
        ENDDO
      ENDIF

      RETURN
      END
C     END OF SUBROUTINE DDMCALCNH3


      SUBROUTINE DDMISRSAVH (H)
      use tracer
      INCLUDE 'camx.prm'
cbk      INCLUDE 'tracer.com'
      INCLUDE 'ddmisrpia.inc'
      REAL*8    H
      SAVH = H
      RETURN
      END
C     END OF SUBROUTINE DDMISRSAVH

