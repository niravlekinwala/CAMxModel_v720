C=====================================================================
C
C   INCLUDE FILE FOR DDM-PM SUBROUTINES FOR RADM
C
C   Dependent upon:
C     tracer.com -> camx.prm  - for MXFDDM
C
C   Log:
C     04/26/2005  bkoo - original development
C     06/25/2008  bkoo - added OMP directive
C     09/19/2013  bkoo - added jTCAC (CaCO3)
C     10/26/2018  bkoo - updated for revised RADM
C
C   Note:
C     Do not modify indices unless you understand this code thoroughly
C
C=====================================================================

      INTEGER    NINPT, NIONS, NSENS
      PARAMETER( NINPT   = 17 )
      PARAMETER( NIONS   = 11 )
      PARAMETER( NSENS   = 24 )

      INTEGER    jTNH3 , jTHNO3, jTFOA , jTCO2 , jTSO2 , jTSVI ,
     &           jTH2O2, jTO3  , jTMHP , jTPAA , jTNACL, jTCAC ,
     &           jTGLY , jTMGLY, jTGLYD, jTHO  , jTSOAC
      PARAMETER( jTNH3   =  1, ! Total (Gas + Dissolved)
     &           jTHNO3  =  2,
     &           jTFOA   =  3,
     &           jTCO2   =  4,
     &           jTSO2   =  5,
     &           jTSVI   =  6,
     &           jTH2O2  =  7,
     &           jTO3    =  8,
     &           jTMHP   =  9,
     &           jTPAA   = 10,
     &           jTNACL  = 11,
     &           jTCAC   = 12,
     &           jTGLY   = 13,
     &           jTMGLY  = 14,
     &           jTGLYD  = 15,
     &           jTHO    = 16,
     &           jTSOAC  = 17 )

      INTEGER    jH     , jNH4   , jNO3   , jHCO2  , jHCO3  , jHSO3  ,
     &           jHSO4  , jOH    , jCO3   , jSO3   , jSO4   , jNH3aq ,
     &           jHNO3aq, jFOAaq , jCO2aq , jSO2aq , jH2O2aq, jO3aq  ,
     &           jMHPaq , jPAAaq , jGLYaq , jMGLYaq, jGLYDaq, jHOaq
      PARAMETER( jH      =  1, ! Ions
     &           jNH4    =  2,
     &           jNO3    =  3,
     &           jHCO2   =  4,
     &           jHCO3   =  5,
     &           jHSO3   =  6,
     &           jHSO4   =  7,
     &           jOH     =  8,
     &           jCO3    =  9,
     &           jSO3    = 10,
     &           jSO4    = 11,
     &           jNH3aq  = 12, ! Undissociated
     &           jHNO3aq = 13,
     &           jFOAaq  = 14,
     &           jCO2aq  = 15,
     &           jSO2aq  = 16,
     &           jH2O2aq = 17,
     &           jO3aq   = 18,
     &           jMHPaq  = 19,
     &           jPAAaq  = 20,
     &           jGLYaq  = 21,
     &           jMGLYaq = 22,
     &           jGLYDaq = 23,
     &           jHOaq   = 24 )

      INTEGER    iKW    , iKA    , iKN    , iKF    , iKC1   , iKC2   ,
     &           iKS1   , iKS2   , iKS3   , iMBNH3 , iMBHNO3, iMBFOA ,
     &           iMBCO2 , iMBSO2 , iMBSVI , iCB
      PARAMETER( iKW     =  1, ! Equilibria
     &           iKA     =  2,
     &           iKN     =  3,
     &           iKF     =  4,
     &           iKC1    =  5,
     &           iKC2    =  6,
     &           iKS1    =  7,
     &           iKS2    =  8,
     &           iKS3    =  9,
     &           iMBNH3  = 10, ! Mass Balances
     &           iMBHNO3 = 11,
     &           iMBFOA  = 12,
     &           iMBCO2  = 13,
     &           iMBSO2  = 14,
     &           iMBSVI  = 15,
     &           iCB     = 16 )! Charge Balance

      REAL       SINI(NINPT,MXFDDM)! Initial Sensitivities of total
      REAL       XMAT(NSENS,MXFDDM)! Solution matrix X
      REAL       Y(NSENS)          ! Current concentrations in the liquid phase
      REAL       SAVSISO2(MXFDDM)  ! Saved SINI(jTSO2)
      REAL       SAVSISVI(MXFDDM)  ! Saved SINI(jTSVI)
      INTEGER    IPAR

      COMMON /DDMRAD/ SINI, XMAT, Y, SAVSISO2, SAVSISVI, IPAR
c$omp threadprivate(/DDMRAD/)

