c----CAMx v7.20 220430
c
c     SOAP.INC contains parameters for OA modules
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        08/25/16 -bkoo-   Revised to include parameters for multiple OA schemes
c        11/28/16 -bkoo-   Revised number of VBS precursors
c        10/16/17 -bkoo-   Added sfacJSOA; removed LSOAP2
c        04/22/19 -gy---   SOAP2.2: Update TRP SOA yields for RO2 autoxidation
c                                   JSOA reduced to 0.001 (Malecha et al. 2018)
c                                   kpolya deactivated because duplicative
c
c-----------------------------------------------------------------------
c     Parameters for SOAP
c-----------------------------------------------------------------------
      integer, parameter :: NPREC_S = 7
      integer, parameter :: NSOAP   = 4   ! No. of CG/SOA pairs
      integer, parameter :: NSOAA   = 2   ! No. of anthropogenic CG/SOA pairs
      integer, parameter :: NNVOL   = 2   ! No. of non-volatile products

      real, parameter :: mwpoa  = 220.    ! MW of POA
      real, parameter :: mwsopa = 220.    ! MW of SOPA
      real, parameter :: mwsopb = 220.    ! MW of SOPB

      real, parameter :: mwsoap(NSOAP) =              ! MW of CG/SOA pairs (g/mol)
     &                                     (/ 150.,   ! CG1
     &                                        150.,   ! CG2
     &                                        180.,   ! CG3
     &                                        180. /) ! CG4
c
c     SOAP1 parameters
c
      real, parameter :: yh_1(3,NPREC_S) = reshape(   ! high-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.315,   0.0186,  0.00679,        ! BNZA
     &              0.0837,  0.0395,  0.0,            ! TOLA
     &              0.0656,  0.0258,  0.00003,        ! XYLA
     &              0.0,     0.282,   0.176,          ! IVOA
     &              0.0791,  0.0132,  0.00115,        ! ISP
     &              0.473,   0.0465,  0.0,            ! TRP
     &              2.14,    0.489,   0.0 /),         ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yl_1(3,NPREC_S) = reshape(   ! low-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.0,     0.0,     0.131,          ! BNZA
     &              0.0,     0.0,     0.125,          ! TOLA
     &              0.0,     0.0,     0.173,          ! XYLA
     &              0.0,     0.317,   0.335,          ! IVOA
     &              0.0791,  0.0132,  0.00115,        ! ISP
     &              0.473,   0.0465,  0.0,            ! TRP
     &              2.14,    0.489,   0.0 /),         ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: csatS_1(NSOAP) =             ! Sat. conc. of CG/SOA pairs (ug/m3)
     &                                     (/ 48.,    ! CG1
     &                                         1.6,   ! CG2
     &                                       140.,    ! CG3
     &                                         2.9 /) ! CG4

      real, parameter :: cstempS_1(NSOAP) =           ! Ref. temp. for sat. conc. (K)
     &                                     (/ 298.,   ! CG1
     &                                        298.,   ! CG2
     &                                        298.,   ! CG3
     &                                        298. /) ! CG4

      real, parameter :: deltahS_1(NSOAP) =           ! Heat of vap. for CG/SOA pairs (J/mol)
     &                                   (/ 20000.,   ! CG1
     &                                      24000.,   ! CG2
     &                                      24000.,   ! CG3
     &                                      57000. /) ! CG4

      real, parameter :: kpolya_1 = 0.035 ! polymerization rate for SOPA [1/hr]; half-life of 20 hr
      real, parameter :: kpolyb_1 = 0.035 ! polymerization rate for SOPB [1/hr]; half-life of 20 hr
      real, parameter :: sfacJSOA_1 = 0.0 ! scaling factor for JSOA: JSOA = sfacJSOA * JNO2
c
c     SOAP2.2 parameters
c
      real, parameter :: yh_2(3,NPREC_S) = reshape(   ! high-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.129,   0.204,   0.0,            ! BNZA
     &              0.187,   0.180,   0.0185,         ! TOLA
     &              0.0593,  0.0343,  0.0123,         ! XYLA
     &              0.0163,  0.319,   0.124,          ! IVOA
     &              0.0286,  0.00850, 0.0,            ! ISP
     &              0.0567,  0.034,   0.0433,         ! TRP
     &              0.105,   0.372,   0.162 /),       ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: yl_2(3,NPREC_S) = reshape(   ! low-NOx yield [ppm/ppm]
c                   More-Vol Less-Vol Non-Vol
     &           (/ 0.253,   0.0870,  0.0,            ! BNZA
     &              0.406,   0.211,   0.109,          ! TOLA
     &              0.206,   0.216,   0.142,          ! XYLA
     &              0.0,     0.388,   0.267,          ! IVOA
     &              0.0588,  0.0108,  0.00335,        ! ISP
     &              0.1133,  0.068,   0.0433,         ! TRP
     &              0.154,   0.453,   0.250 /),       ! SQT
     &           (/ 3, NPREC_S /) )

      real, parameter :: csatS_2(NSOAP) =             ! Sat. conc. of CG/SOA pairs (ug/m3)
     &                                    (/ 14.,     ! CG1
     &                                        0.31,   ! CG2
     &                                       26.,     ! CG3
     &                                        0.45 /) ! CG4

      real, parameter :: cstempS_2(NSOAP) =           ! Ref. temp. for sat. conc. (K)
     &                                     (/ 300.,   ! CG1
     &                                        300.,   ! CG2
     &                                        300.,   ! CG3
     &                                        300. /) ! CG4

      real, parameter :: deltahS_2(NSOAP) =           ! Heat of vap. for CG/SOA pairs (J/mol)
     &                                  (/ 116000.,   ! CG1
     &                                     147000.,   ! CG2
     &                                     118000.,   ! CG3
     &                                     123000. /) ! CG4

      real, parameter :: kpolya_2 = 0.0   ! polymerization rate for SOPA [1/hr]; deactivated
      real, parameter :: kpolyb_2 = 0.0   ! polymerization rate for SOPB [1/hr]; deactivated
      real, parameter :: sfacJSOA_2 = 0.001 ! scaling factor for JSOA: JSOA = sfacJSOA * JNO2

      real :: kpolya, kpolyb, sfacJSOA
      real, dimension(NSOAP) :: csatS, cstempS, deltahS
      real, dimension(3,NPREC_S) :: y_h, y_l

      common /soapcom/ kpolya, kpolyb, sfacJSOA, csatS, cstempS, deltahS,
     &                 y_h, y_l

c-----------------------------------------------------------------------
c     Parameters for VBS
c-----------------------------------------------------------------------
c
      integer, parameter :: NPREC_V = 10
      integer, parameter :: NVBS    = 4   ! No. of volatility bins (volatile bins only)
      integer, parameter :: NSET    = 5   ! No. of basis sets

      real, parameter :: mwPAS0 = 172.    ! MW of PAS0
      real, parameter :: mwPBS0 = 172.    ! MW of PBS0
      real, parameter :: mwPAP0 = 278.    ! MW of PAP0
      real, parameter :: mwPCP0 = 278.    ! MW of PCP0
      real, parameter :: mwPFP0 = 205.    ! MW of PFP0

      real, parameter :: mwVBS(NVBS*NSET) = 
     &                                   (/ 167., 163., 158., 153.,   ! PAS1-PAS4
     &                                      167., 163., 158., 153.,   ! PBS1-PBS4
     &                                      275., 272., 268., 266.,   ! PAP1-PAP4
     &                                      275., 272., 268., 266.,   ! PCP1-PCP4
     &                                      208., 211., 213., 215. /) ! PFP1-PFP4

      real, parameter :: csatV(NVBS*NSET) =
     &                                   (/ 1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000.,
     &                                      1., 10., 100., 1000. /)

      real, parameter :: cstempV(NVBS*NSET) =
     &                                   (/ 298., 298., 298., 298.,
     &                                      298., 298., 298., 298.,
     &                                      298., 298., 298., 298.,
     &                                      298., 298., 298., 298.,
     &                                      298., 298., 298., 298. /)

      real, parameter :: deltahV(NVBS*NSET) =
     &                             (/ 35000., 35000., 35000., 35000.,
     &                                35000., 35000., 35000., 35000.,
     &                                85000., 74000., 63000., 52000.,
     &                                85000., 74000., 63000., 52000.,
     &                                85000., 81000., 77000., 73000. /)

c-----------------------------------------------------------------------
c     Parameters for common arrays
c-----------------------------------------------------------------------
      integer, parameter :: NOXID = 6                     ! No. of oxidant used in precursor chemistry
      integer, parameter :: NPREC = MAX(NPREC_S, NPREC_V) ! No. of SOA precursors
      integer, parameter :: NCG = MAX(NSOAP, NVBS*NSET)   ! No. of CG species
      integer, parameter :: NNV = MAX(NNVOL, NSET)        ! No. of non-volatile products

