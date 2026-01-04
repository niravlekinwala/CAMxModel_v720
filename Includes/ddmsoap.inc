C=====================================================================
C
C   INCLUDE FILE FOR DDM-PM SUBROUTINES FOR SOAP
C
C   Dependent upon:
C     tracer.com -> camx.prm  - for MXFDDM
C
C   Log:
C     05/20/2005  bkoo - original development
C     06/25/2008  bkoo - revised for CAMx v4.5 SOA module update
C                      - added OMP directive
C     08/25/2016  bkoo - updated NSENS for new SOAP
C
C=====================================================================

      INTEGER    NSENS
      PARAMETER( NSENS   =  4 )    ! MUST BE CONSISTENT WITH NSOAP IN SOAP.INC

      REAL       SINI(NSENS,MXFDDM)! Initial Sensitivities of total
      REAL       SPRE(MXFDDM)      ! Initial Sensitivities of pre-existing organics in mole
      REAL       SIPA(MXFDDM)      ! Initial Sensitivities of polymerized anthro SOA
      REAL       SIPB(MXFDDM)      ! Initial Sensitivities of polymerized bio SOA
      REAL       XMAT(NSENS,MXFDDM)! Solution matrix X
      INTEGER    IPAR

      COMMON /DDMSOA/ SINI, SPRE, SIPA, SIPB, XMAT, IPAR
c$omp threadprivate(/DDMSOA/)

