c----CAMx v7.20 220430
c
c     VBS.INC contains common parameters and variables for VBS
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        12/07/14 -bkoo-      Created
c        01/08/16 -bkoo-      Updated for Revised VBS
c        08/25/16 -bkoo-      Modified to contain emiss-related params only
c
c-----------------------------------------------------------------------
c
      logical, parameter :: LVBSPREPROC = .TRUE.
c
      integer, parameter :: NVOLBIN = 4   ! MUST BE CONSISTENT WITH NVBS IN SOAP.INC
c
c     Volatility distributions for sector-specific POA emissions
c
      real, parameter, dimension(0:NVOLBIN) ::
     &                      POA_GV_EF = (/0.27,0.15,0.26,0.15,0.17/),
     &                      POA_DV_EF = (/0.03,0.25,0.37,0.24,0.11/),
     &                      POA_OP_EF = (/0.09,0.09,0.14,0.18,0.50/),
     &                      POA_BB_EF = (/0.2, 0.1, 0.1, 0.2, 0.4 /),
     &                      POA_MC_EF = (/0.35,0.35,0.10,0.10,0.10/)
c
c     Common Variables
c
      logical :: lvbs
      integer, dimension(0:NVOLBIN) :: kpap_c, kpcp_c, kpfp_c

      common /vbscom/ lvbs, kpap_c, kpcp_c, kpfp_c

