C=====================================================================
C
C   INCLUDE FILE FOR DDM-PM SUBROUTINES FOR ISORROPIA
C
C   Dependent upon:
C     tracer.com -> camx.prm  - for MXFDDM
C
C   Log:
C     03/03/2005  bkoo - original development
C     06/25/2008  bkoo - added OMP directive
C
C   Note:
C     In general, index parameters are consistent with those of ISORROPIA
C     Do not modify them unless you understand this code thoroughly
C
C=====================================================================

      INTEGER    NINPT, NIONSP1, NSENS
      PARAMETER( NINPT   =  5 )! = NCOMP     (NCOMP -> ISRPIA.INC)
      PARAMETER( NIONSP1 =  8 )! = NIONS + 1 (NIONS -> ISRPIA.INC)
      PARAMETER( NSENS   = 20 )! Number of the sensitivity coefficients

      INTEGER    jTNA   , jTSO4  , jTNH4  , jTNO3  , jTCL
      PARAMETER( jTNA    =  1, ! Total (Gas + Aerosol)
     &           jTSO4   =  2,
     &           jTNH4   =  3,
     &           jTNO3   =  4,
     &           jTCL    =  5 )

      INTEGER    jH     , jNA    , jNH4   , jCL    , jSO4   , jHSO4  ,
     &           jNO3   , jH2O   , jNH3   , jHCL   , jHNO3  , jNH42S4,
     &           jNH4HS4, jNACL  , jNA2SO4, jNANO3 , jNH4NO3, jNH4CL ,
     &           jNAHSO4, jLC
      PARAMETER( jH      =  1, ! IONS
     &           jNA     =  2,
     &           jNH4    =  3,
     &           jCL     =  4,
     &           jSO4    =  5,
     &           jHSO4   =  6,
     &           jNO3    =  7,
     &           jH2O    =  8, ! WATER
     &           jNH3    =  9, ! GASES
     &           jHCL    = 10,
     &           jHNO3   = 11,
     &           jNH42S4 = 12, ! SOLIDS
     &           jNH4HS4 = 13,
     &           jNACL   = 14,
     &           jNA2SO4 = 15,
     &           jNANO3  = 16,
     &           jNH4NO3 = 17,
     &           jNH4CL  = 18,
     &           jNAHSO4 = 19,
     &           jLC     = 20 )

      INTEGER    mNACL  , mNA2SO4, mNANO3 , mNH42S4, mNH4NO3, mNH4CL ,
     &           mH2SO4 , mHHSO4 , mNH4HS4, mHNO3  , mHCL   , mNAHSO4,
     &           mLC
      PARAMETER( mNACL   =  1, ! ION PAIRS
     &           mNA2SO4 =  2,
     &           mNANO3  =  3,
     &           mNH42S4 =  4,
     &           mNH4NO3 =  5,
     &           mNH4CL  =  6,
     &           mH2SO4  =  7,
     &           mHHSO4  =  8,
     &           mNH4HS4 =  9,
     &           mHNO3   = 10,
     &           mHCL    = 11,
     &           mNAHSO4 = 12,
     &           mLC     = 13 )

      INTEGER    iK1    , iK2    , iK3    , iK4    , iK5    , iK6    ,
     &           iK7    , iK8    , iK9    , iK10   , iK11   , iK12   ,
     &           iK13   , iMBNA  , iMBSO4 , iMBNH4 , iMBNO3 , iMBCL  ,
     &           iCB    , iZSR
      PARAMETER( iK1     =  1, ! Equilibrium reactions
     &           iK2     =  2,
     &           iK3     =  3,
     &           iK4     =  4,
     &           iK5     =  5,
     &           iK6     =  6,
     &           iK7     =  7,
     &           iK8     =  8,
     &           iK9     =  9,
     &           iK10    = 10,
     &           iK11    = 11,
     &           iK12    = 12,
     &           iK13    = 13,
     &           iMBNA   = 14, ! Mass Balances
     &           iMBSO4  = 15,
     &           iMBNH4  = 16,
     &           iMBNO3  = 17,
     &           iMBCL   = 18,
     &           iCB     = 19, ! Charge Balance
     &           iZSR    = 20 )! ZSR relationship

      REAL       SINI(NINPT,MXFDDM)! Initial Sensitivities of total
      REAL       XMAT(NSENS,MXFDDM)! Solution matrix X
      REAL       DALFA(3,MXFDDM)   ! HCL; HNO3; NH3
      REAL*8     SAVH              ! H+ before CALCNH3

      LOGICAL    LMDRH         ! Logical flag for MDRH region
      INTEGER    IPAR

      COMMON /DDMISR/ SINI, XMAT, DALFA, SAVH, LMDRH, IPAR
c$omp threadprivate(/DDMISR/)

